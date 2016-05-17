#include "./thread.h"

#include <alloca.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iostream>

using namespace std;

using namespace Threads;

// By default, threads block all signals. This flag disables signal blocking.
static bool disable_thread_signal_block = ::getenv("MT_DISABLE_THREAD_SIGNAL_BLOCK");

// This code is used to determine the size of the TLS space.
// TLS variables live in the stack space allocated for each thread, so TLS size
// is used to adjust thread stack size.
// This only works on GCC with GLIBC because it relies on an internal library
// symbol _dl_get_tls_static_info.
// TODO: This implementation does not work for dynamically loaded libraries,
//       where we need to call dlsym(TLD_NEXT, "_dl_get_tls_static_info") to
//       get the address if this symbol.
#ifdef __GLIBC__
#ifdef __i386__
#define GET_TLS_ATTR __attribute__((weak)) __attribute ((regparm (3), stdcall))
#else
#define GET_TLS_ATTR __attribute__((weak))
#endif
extern "C" void _dl_get_tls_static_info(size_t*, size_t*) GET_TLS_ATTR;
#undef GET_TLS_ATTR
static size_t GetTLSSize() {
    size_t tls_size = 0;
    size_t tls_align;
    _dl_get_tls_static_info(&tls_size, &tls_align);
    return tls_size;
}
#else
static size_t GetTLSSize() {
    return 0;
}
#endif // __GLIBC__
size_t xxx = GetTLSSize();

// Default system stack size in bytes. This value should be empirically adjusted.
// TODO: this should be thread-safe!!!
//static AtomicGuard<size_t> default_stack_size(2048*1024);
static size_t default_stack_size(2048*1024);

// Some platforms do not define PTHREAD_STACK_MIN
#ifndef PTHREAD_STACK_MIN
#define PTHREAD_STACK_MIN 16384
#endif // PTHREAD_STACK_MIN 

// Return the smallest value divisible by alignment that is greater than the
// specified address.
static inline uintptr_t align_up(uintptr_t address, size_t alignment) {
    uintptr_t mask = alignment - 1;
    return (address + mask) & ~mask;
}

// Return the minimum valid stack size that can be passed to
// pthread_attr_setstacksize().
static size_t MinStackSize(size_t stack_size) {
  if (stack_size < PTHREAD_STACK_MIN) stack_size = PTHREAD_STACK_MIN;

  // Make stack size a a multiple of the system page size.
  static size_t page_size = sysconf(_SC_PAGESIZE);
  return align_up(stack_size, page_size);
}

void Thread::Options::set_default_stack_size(size_t stack_size) {
    stack_size = MinStackSize(stack_size);
    ::default_stack_size = stack_size;
}

size_t Thread::Options::default_stack_size() {
    return ::default_stack_size;
}

// Default system stack offset in bytes. This value should be empirically adjusted.
// The goal is to avoid L1-cache and TLB conflicts, the value of 1K per thread
// plus 256 bytes is recommended. The default value is good for up to 16
// threads.
// TODO: this should be thread-safe!!!
//static AtomicGuard<size_t> default_stack_offset(16*1024 + 256);
static size_t default_stack_offset(16*1024 + 256);

void Thread::Options::set_default_stack_offset(size_t stack_offset) {
    ::default_stack_offset = stack_offset;
}

size_t Thread::Options::default_stack_offset() {
    return ::default_stack_offset;
}

// Default guard size is not configurable, just hard-wired to some empirical values:
// 1MB guard region on 64-bit systems and 16kB guard region on 32-bit systems.
static size_t default_guard_size() {
    return (sizeof(void*) > 4) ? (1 << 20) : (1 << 14);
}

// Thread options.
Thread::Options::Options()
    : stack_size_(0),
      stack_offset_(0),
      guard_size_(0),
      joinable_(true)
{
}


Thread::Thread(const Thread::Options& options)
    : options_(options),
      tid_(0),
      created_(false),
      needs_join_(false)
{
}

Thread::~Thread() {
    if (needs_join_) {
        LOG_ERROR_NONFATAL << "Joinable thread was not joined - memory and other resources will be leaked!";
    }
}

// Create the thread and run the payload.
void Thread::Start() {
    CHECK(!created_) << "Start() called on a running thread!";
    created_ = true;
    needs_join_ = options_.joinable();

    // Prepare thread attributes.
    pthread_attr_t attr;
    CHECK_EQ(pthread_attr_init(&attr), 0);
    int detach = options_.joinable() ? PTHREAD_CREATE_JOINABLE : PTHREAD_CREATE_DETACHED;
    CHECK_EQ(pthread_attr_setdetachstate(&attr, detach), 0);

    // Set up stack size. 
    // Use the system default size if stack size is not specified.
    size_t stack_size = Thread::Options::default_stack_size();
    if (options_.stack_size() != 0) {
        stack_size = options_.stack_size() + options_.stack_offset();
    }

    // Set up guard size.
    // Use the system default size if guard size is not specified.
    size_t guard_size = ::default_guard_size();
    if (options_.guard_size() != 0) {
        guard_size = options_.guard_size();
    }

  // GLIBC (at least one that use NPTL) takes the guard size out of the space
  // allocated for the stack, rather than adding the guard on to the end.  Make
  // up for that here.
#ifdef __GLIBC__
    stack_size += guard_size;
#endif // __GLIBC__

    // Padd the stack by the TLS size, make sure the value is valid, and set the size.
    stack_size = MinStackSize(stack_size + GetTLSSize());
    CHECK_EQ(pthread_attr_setstacksize(&attr, stack_size), 0)
        << ": specified stack size = " << stack_size
        << ", PTHREAD_STACK_MIN= " << PTHREAD_STACK_MIN;

    // Set the guard size.
    CHECK_EQ(pthread_attr_setguardsize(&attr, guard_size), 0);

    // Time to do the real work!
    int rc = pthread_create(&tid_, &attr, Payload, this);
    CHECK_EQ(rc, 0) << "Failed to start a thread: " << strerror(rc)
        << ((rc == ENOMEM) ? " Thread stack size may be too large." : "" );

    // Clean up temporary objects.
    CHECK_EQ(pthread_attr_destroy(&attr), 0);
}

// Join the thread (if it is joinable).
void Thread::Join() {
    CHECK(options_.joinable()) << "Attempting to join a detached thread!";
    CHECK(created_) << "Attempting to join a non-existing thread!";
    if (!needs_join_) return;
    int rc = pthread_join(tid_, NULL);
    CHECK_EQ(rc, 0) << "Failed to join a thread: " << strerror(rc)
        << ((rc == EDEADLK) ? " Deadlock, is the thread joining itself?" : "");
    needs_join_ = false;
}

// Cancel the thread.
void Thread::Cancel() {
    CHECK(created_) << "Attempting to cancel a non-existing thread!";
    int rc = pthread_cancel(tid_);
    CHECK_EQ(rc, 0) << "Failed to cancel a thread: " << strerror(rc);
}

// Thread payload, this is where actual work is done.
void* Thread::Payload(void* arg) {
    // Make sure that this thread does not get any signals.
    sigset_t set;
    if (!disable_thread_signal_block) { // Default, most signals are disabled
        sigfillset( &set );
        sigdelset( &set, SIGPROF );
        sigdelset( &set, SIGSEGV );
        sigdelset( &set, SIGBUS );
    } else { // Make sure that this thread does not get the SIGCLD signal from MGLS
        sigemptyset( &set );
        sigaddset( &set, SIGCHLD );
    }
    pthread_sigmask(SIG_BLOCK, &set, NULL);

    // Recover thread pointer.
    Thread* this_thread = static_cast<Thread*>(arg);

    // Stack size in Start() was padded by the specified amount. Now reserve that
    // memory to avoid cache and TSL conflicts.
    const char* reserved = NULL;
    if (this_thread->options_.stack_offset() > 0) {
        reserved = reinterpret_cast<char*>(alloca(this_thread->options_.stack_offset()));
    }

    // And now, the moment we have all been waiting for...
    this_thread->Run();

    // The Thread object may have been deleted by the Run() method, so it
    // should not be accessed beyond this point.
    this_thread = NULL;

    return NULL;
}


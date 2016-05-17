
all: DEFAULT
thread_lib =  ./pthread_library
PROG = main
thread_o = $(thread_lib)/common.o $(thread_lib)/thread.o
DEFAULT: $(PROG)


CCPATH   = /usr/bin/g++
 


# Compile the application code
objects = $(PROG).o 
CCFLAG = -g -O3 #-pg

$(PROG).o: $(PROG).cpp Mythread.h
	$(CCPATH) $(CCFLAG) -o $(PROG).o \
	 -c $(PROG).cpp 


Mythread.o : Mythread.h Mythread.cpp Demagnetization_factors.cpp $(thread_lib)/thread.h
	$(CCPATH) $(CCFLAG) -o Mythread.o \
	-c Mythread.cpp
	

# Link the executable

$(PROG): $(PROG).o $(thread_o)  Mythread.o
	$(CCPATH) $(CCFLAG)  -lpthread -o $(PROG)  $(PROG).o $(thread_o)  Mythread.o
          





clean: 
	@/bin/rm  -rf $(PROG).o 
	@/bin/rm  -rf $(PROG)
	@/bin/rm -rf *.o
	cd $(thread_lib) && $(MAKE) clean

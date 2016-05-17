# Except as specified in the OpenAccess terms of use of Cadence or Silicon
# Integration Initiative, this material may not be copied, modified,
# re-published, uploaded, executed, or distributed in any way, in any medium,
# in whole or in part, without prior written permission from Cadence.
#
#                Copyright 2002-2005 Cadence Design Systems, Inc.
#                           All Rights Reserved.
#
#  $Author: pdsim $
#  $Revision: #21 $
#  $Date: 2006/11/15 $
# ******************************************************************************
# ******************************************************************************


#include ../macro.defs

abspath = /app/design/puneet/projects/Shaodiwang/MeRAM/MTJ_variation/Precessional/Debug/Pthreads_version/pthread_library
objects = $(abspath)/common.o $(abspath)/thread.o
DEFAULT: $(objects)



CCPATH   = /usr/bin/g++

CCFLAG = #-O4 #-pg

$(objects): %.o:%.C
	$(CCPATH) $(CCFLAG) -o $@ \
	-c $<


	

	

# Link the executable

#$(PROG): $(PROG).o common.o timers.o thread.o Mythread.o
#	$(CCPATH) $(CCFLAG)  -lrt -lpthread -o $(PROG)  $(PROG).o Mythread.o thread.o timers.o common.o
          



#$(PROG): $(PROG).o  $(OA_LIB_LIST)
#	$(CCPATH) $(DEBUG) $(CXXOPTS) -o $(PROG)  $(PROG).o \
#	 $(CCLNFLAG) 


clean: 
	@/bin/rm  -rf $(PROG).o 
	@/bin/rm  -rf $(PROG)
	@/bin/rm -rf *.o

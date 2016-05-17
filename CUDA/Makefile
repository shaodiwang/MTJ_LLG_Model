
all : DEFAULT
PROG = WERSim
NVCC = /u/local/cuda/5.0/bin/nvcc
DEFAULT: $(PROG)
CUPROG = LLG 

CCFLAG = -I/u/local/cuda/5.0/include/ -arch=sm_13 --optimize 2 #-G -g#/u/local/cuda/current/include/ #-G -g
CCPATH   = g++

# Compile the application code
objects = $(PROG).o 


$(PROG).o: $(PROG).cu
	$(NVCC) $(CCFLAG) -o $(PROG).o \
	 -c $(PROG).cu




# Link the executable

$(PROG): $(PROG).o 
	$(NVCC) $(CCFLAG)  -o $(PROG)  $(PROG).o
          



#$(PROG): $(PROG).o  $(OA_LIB_LIST)
#	$(CCPATH) $(DEBUG) $(CXXOPTS) -o $(PROG)  $(PROG).o \
#	 $(CCLNFLAG) 


clean: 
	@/bin/rm  -rf $(PROG).o 
	@/bin/rm  -rf $(PROG)
	@/bin/rm  -rf $(CUPROG).o

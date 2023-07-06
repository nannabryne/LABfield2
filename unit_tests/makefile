# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# NB! Not to be used in this directory. 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# makefile for running simple unit tests (LATfield2)


# programming environment
COMPILER     := mpigxx
INCLUDE      := -I../../../LATfield2/# path to LATfield2 
LIB          := -lfftw3_omp -lfftw3 -lm -lhdf5
COMPILERFLAGS:= -fopenmp

# targets and sources
EXEC         := main
SOURCE       := main.cpp 
HEADERS      := ../utils.h
TARGETS 	 := main



# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5#-DH5_HAVE_PARALLEL


# optional compiler settings 
DPROGRAM 	 :=






# further compiler options
OPT          := -O3 -std=c++11




# parallel grid:
n?=2
m?=2
np=$(shell echo ${n}*${m} | bc)



bench_this?=n# whether to make the original sample (y) or do the analysis (n)
ifeq ($(strip $(bench_this)), y)
DPROGRAM +=-D_BENCH
n=4
m=10
endif



# 	-	-	-	-	-	-	-	-	-	
#	COMMANDS
# 	-	-	-	-	-	-	-	-	-


# 	[make]
$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $(COMPILERFLAGS) $< -o $@ $(OPT) $(DLATFIELD2) $(DPROGRAM) $(INCLUDE) $(LIB)


# command to build C++ program:
all: 
	$(COMPILER) $(COMPILERFLAGS) $(SOURCE) -o $(EXEC) $(OPT) $(DLATFIELD2) $(DPROGRAM) $(INCLUDE) $(LIB)
	mpirun -np 4 -ppn 1 ./$(EXEC) -n 2 -m 2


run: 
	mpirun -np $(np) ./$(EXEC) -n $(n) -m $(m)


bench: $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) -D_BENCH $(INCLUDE) $(LIB)
	mpirun -np 80 -ppn 20 ./$@ -n 4 -m 20
	rm $@



# command to delete object files and executables
.PHONY: clean
clean:
	rm -f $(wildcard *.o) $(EXEC) 




# testing alternative approach


test: 
	cd fft_execution && make && make run
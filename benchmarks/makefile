# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# NB! Not to be used in this directory. 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# makefile for running simple unit tests (LATfield2)


# programming environment
COMPILER     := mpigxx
INCLUDE      := -I../../../LATfield2/# path to LATfield2 
LIB          := -lfftw3 -lm -lhdf5


# targets and sources
EXEC         := main
SOURCE       := main.cpp 
HEADERS      := ../utils.h
TARGETS 	 := main


# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5 #-DH5_HAVE_PARALLEL


# optional compiler settings (gevolution)
DPROGRAM 	 :=


# further compiler options
OPT          := -O3 -std=c++11




# parallel grid:
n?=2
m?=2
np=$(shell echo ${n}*${m} | bc)


# 	-	-	-	-	-	-	-	-	-	
#	COMMANDS
# 	-	-	-	-	-	-	-	-	-


# 	[make]
$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DPROGRAM) $(INCLUDE) $(LIB)


# command to build C++ program:
all: $(EXEC) 


run: 
	mpirun -np $(np) ./$(EXEC) -n $(n) -m $(m)



# command to delete object files and executables
.PHONY: clean
clean:
	rm -f $(wildcard *.o) $(EXEC)




# testing alternative approach


test: 
	cd fft_execution && make && make run
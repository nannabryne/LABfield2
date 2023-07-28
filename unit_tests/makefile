

# programming environment
COMPILER  		:= mpigxx
INCLUDE     	:= -I../../LATfield2/# path to LATfield2
LIB         	:= -lfftw3_omp -lfftw3 -lm -lhdf5
COMPILERFLAGS	:= -fopenmp


# targets and sources
EXEC 			:= main.exe
SOURCE			:= main.cpp


# flags for LATfield2
DLATFIELD2  	:= -DFFT3D -DHDF5


# flags for this program
DPROGRAM 		:= -DLOOPCORR -DPROJECTION -DFOURIER -DPART


# further compiler options
OPT         	:= -O3 -std=c++11






# parallel grid:
n?=4
m?=4
o?=20
np=$(shell echo ${n}*${m} | bc)


ifeq ($(strip $(o)), 1)
COMPILERFLAGS =
LIB = -lfftw3 -lm -lhdf5
endif



$(EXEC): $(SOURCE) makefile
	$(COMPILER) $(COMPILERFLAGS) $< -o $@ $(OPT) $(DLATFIELD2) $(DPROGRAM) $(INCLUDE) $(LIB)



run: 
	mpirun -np $(np) ./$(EXEC) -n $(n) -m $(m)



.PHONY: clean
clean:
	rm -f $(EXEC) $(wildcard *.o) 
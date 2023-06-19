
# makefile for running simple unit tests (LATfield2)

# programming environment
COMPILER     := mpigxx
INCLUDE      := -I../../LATfield2/# path to LATfield2 
LIB          := -lfftw3 -lm -lhdf5


# targets and sources
SOURCE_DIR   := unit_tests
BUILD_DIR	 := build
INCLUDE_DIR  := boilerplate_examples

EXEC         := $(BUILD_DIR)/program
SOURCE       := $(SOURCE_DIR)/field_manipulation.cpp#for now
HEADERS      := $(INCLUDE_DIR)/field_manipulation.hpp
TARGETS 	 := program


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





# command to create build directory:
dir:
	mkdir -p $(BUILD_DIR)

# command to build all C++ programs:
all: $(EXEC) 



# mpigxx -O3 -std=c++11 -o $(EXEC) $(SOURCE) $(INCLUDE) -DFFT3D -DHDF5 -lfftw3 -lm -lhdf5


run: 
	mpirun -np $(np) $(EXEC) -n $(n) -m $(m)



# command to delete object files and executables
.PHONY: clean
clean:
	rm -f $(BUILD_DIR)/$(wildcard *.o) $(EXEC)
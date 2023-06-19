
# makefile for running simple unit tests (LATfield2)

# programming environment
COMPILER     := mpigxx
INCLUDE      := -I../../../LATfield2/# path to LATfield2 
LIB          := -lfftw3 -lm -lhdf5


# targets and sources
SOURCE_DIR   := boilerplate_examples/ 
BUILD_DIR	 := build/

EXEC         := program
SOURCE       := $(SOURCE_DIR)main.cpp
HEADERS      := $(wildcard *.hpp)


# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5


# optional compiler settings (gevolution)
DPROGRAM 	 :=


# further compiler options
OPT          := -O3 -std=c++11

# 	[make]
$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DPROGRAM) $(INCLUDE) $(LIB)




# parallel grid:
n?=2
m?=2
np=$(shell echo ${n}*${m} | bc)



# 	-	-	-	-	-	-	-	-	-	
#	COMMANDS
# 	-	-	-	-	-	-	-	-	-

# command to create build directory:
dir:
	mkdir -p $(BUILD_DIR)

# command to build all C++ programs:
all: $(BUILD_DIR)/$(TARGETS)



# mpigxx -O3 -std=c++11 -o $(EXEC) $(SOURCE) $(INCLUDE) -DFFT3D -DHDF5 -lfftw3 -lm -lhdf5


run: 
	mpirun -np $(np) ./$(EXEC) -n $(n) -m $(m)



clean:
	rm $(EXEC) $(wildcard test_*)
# A simple Makefile for CDFCI

# compiler path (icc, icpc, g++) (icc highly recommended)
CC = g++
# compiler flag
CCFLAG = -O3 -std=c++14
# compiler debug flag
DBGFLAG = -g -Wall
# compiler OpenMP flag
OMP_FLAG = -fopenmp
# Mac OS X uses -openmp
ifneq ($(OS), Windows_NT)
ifeq ($(shell uname -s), Darwin)
        OMP_FLAG = -openmp
endif
endif
# Intel compiler uses -qopenmp
ifeq (icc, $(findstring icc, $(CC)))
        OMP_FLAG = -qopenmp
endif
ifeq (icpc, $(findstring icpc, $(CC)))
        OMP_FLAG = -qopenmp
endif

# executable name
# OpenMP disabled, single thread
TARGET = cdfci
# OpenMP enabled
TARGET_OMP = $(TARGET)_omp

# source file containing the main() function
MAIN_SOURCE = src/main.cpp

# build test
# compiler
# Use g++ for CI integration
CC_TEST = $(CC)
OMP_TEST_FLAG = $(OMP_FLAG)
# test directory
DIR_TEST = test
MAIN_TEST = src/test_system.cpp
TARGET_TEST = cdfci_test
TARGET_OMP_TEST = $(TARGET_TEST)_omp

# default rule
default: all

.PHONY: all
all:
	$(CC) $(CCFLAG) $(MAIN_SOURCE) -o $(TARGET)
	$(CC) $(CCFLAG) $(OMP_FLAG) $(MAIN_SOURCE) -o $(TARGET_OMP)

.PHONY: debug
debug:
	$(CC) $(CCFLAG) $(DBGFLAG) $(MAIN_SOURCE) -o $(TARGET)
	$(CC) $(CCFLAG) $(DBGFLAG) $(OMP_FLAG) $(MAIN_SOURCE) -o $(TARGET_OMP)

.PHONY: build_test
build_test:
	$(CC_TEST) $(CCFLAG) $(DIR_TEST)/$(MAIN_TEST) -o $(DIR_TEST)/$(TARGET_TEST)
	$(CC_TEST) $(CCFLAG) $(OMP_TEST_FLAG) $(DIR_TEST)/$(MAIN_TEST) -o $(DIR_TEST)/$(TARGET_OMP_TEST)

.PHONY: test
test: build_test
	cd $(DIR_TEST) && ./$(TARGET_TEST) -D
ifeq ($(OS), Windows_NT)
	set OMP_NUM_THREADS=2
else
	export OMP_NUM_THREADS=2
endif
	cd $(DIR_TEST) && ./$(TARGET_OMP_TEST) -D

.PHONY: clean
clean:
	rm -f $(TARGET)
	rm -f $(TARGET_OMP)
	rm -f $(DIR_TEST)/$(TARGET_TEST)
	rm -f $(DIR_TEST)/$(TARGET_OMP_TEST)

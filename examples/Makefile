# Compiler options
CXX ?= g++
CXXFLAGS := -O3 -ffast-math -march=native -fopenmp
LDFLAGS = -L../src
LDLIBS = -lLBM
INCLUDE = -I../include

# Source 
LOC = bin
SRC := naive_implementation.cpp naive_improved_implementation.cpp naive_parallel.cpp
PROGS := $(LOC)/naive_implementation $(LOC)/naive_improved_implementation $(LOC)/naive_parallel

all: $(PROGS)

$(LOC)/naive_implementation: naive_implementation.cpp
	mkdir -p $(LOC)
	$(CXX) $(INCLUDE) $(LDFLAGS) $(LDLIBS) $(CXXFLAGS) naive_implementation.cpp -o $@


$(LOC)/naive_improved_implementation: naive_improved_implementation.cpp
	mkdir -p $(LOC)
	$(CXX) $(INCLUDE) $(LDFLAGS) $(LDLIBS) $(CXXFLAGS) naive_improved_implementation.cpp -o $@

$(LOC)/naive_parallel: naive_parallel.cpp
	mkdir -p $(LOC)
	$(CXX) $(INCLUDE) $(LDFLAGS) $(LDLIBS) $(CXXFLAGS) naive_parallel.cpp -o $@


clean:
	rm -rf $(LOC)




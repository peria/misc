CXX=nvcc
CXXFLAGS=-O2

.PHONY: all
all: bbp

bbp: bbp.cu
	$(CXX) $(CXXFLAGS) -o $@ $^


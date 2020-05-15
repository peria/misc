CXX=clang++
CXXFLAGS=-O2 -std=c++17

.PHONY: all
all: bbp

bbp: bbp.cu bbp.cc
	$(CXX) $(CXXFLAGS) -o $@ bbp.cc


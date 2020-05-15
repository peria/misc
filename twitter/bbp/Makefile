CXX=clang++
CXXFLAGS=-O2 -std=c++17 -W

.PHONY: all
all: bbp

bbp: bbp.cu bbp.cc Makefile
	$(CXX) $(CXXFLAGS) -o $@ bbp.cc


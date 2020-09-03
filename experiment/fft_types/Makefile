CXX=clang++
CXXFLAGS=-O2 -std=c++14

HEADS=stockham_dit.h stockham_dif.h cooley_ditf.h cooley_dift.h
HEADS+=fft.h complex.h

all: benchmark

benchmark: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@
	./$@

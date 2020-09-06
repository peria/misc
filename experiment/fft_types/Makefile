CXX=clang++
CXXFLAGS=-O2 -std=c++14

HEADS=stockham_dit.h stockham_dif.h cooley_dit1.h cooley_dif1.h cooley_dit2.h cooley_dif2.h
HEADS+=fft.h complex.h

all: benchmark

benchmark: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@
	./$@

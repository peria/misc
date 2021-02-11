CXX=clang++
CXXFLAGS=-O2 -std=c++14

#CXXFLAGS+=-O0 -g

HEADS=fft.h complex.h
HEADS+=pmp.h
HEADS+=cooley_dif.h
#stockham_dit.h stockham_dif.h cooley_dit1.h cooley_dif1.h cooley_dit2.h cooley_dif2.h

all: run

run: benchmark
	./benchmark

benchmark: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@

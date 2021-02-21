CXX=clang++
CXXFLAGS=-O2 -std=c++14

#CXXFLAGS+=-O0 -g

HEADS=fft.h complex.h
HEADS+=pmp.h
HEADS+=pmp2.h
HEADS+=pmp5.h
HEADS+=cooley.h
HEADS+=stockham_dit.h
HEADS+=stockham_dif.h
HEADS+=stockham6.h

all: run

run: benchmark
	./benchmark

benchmark: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@

CXX=clang++
CXXFLAGS=-O2 -std=c++14

#CXXFLAGS+=-O0 -g

HEADS=fft.h complex.h
HEADS+=pmp.h
HEADS+=dpmp.h
HEADS+=pmp2.h
HEADS+=pmp5.h
HEADS+=cooley.h
HEADS+=stockham_dit.h
HEADS+=stockham_dif.h
HEADS+=stockham6.h

all: run

run: rft_bench
	./rft_bench

rft_bench: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) -DREAL_FFT=1 $< -o $@

dft_bench: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@

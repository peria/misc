CXX=clang++
CXXFLAGS=-O2 -std=c++11 -march=native

all: fmt

fmt: fft.cc ntt.cc mtt.cc main.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

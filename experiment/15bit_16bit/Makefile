CXX=clang++
CXXFLAGS=-O2 -std=c++14 -W

all: test

test: main.cc fft.cc fft.h base.h
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.cc, $^)

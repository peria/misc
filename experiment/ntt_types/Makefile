CXX=clang++
CXXFLAGS=-O2 -std=c++14

#CXXFLAGS+=-O0 -g

HEADS=$(shell ls *.h)
.PHONY: all run check

all: run

run: ntt_bench
	./ntt_bench

check: ntt_bench
	./ntt_bench check

ntt_bench: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@

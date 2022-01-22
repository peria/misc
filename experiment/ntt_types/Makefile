CXX=clang++
CXXFLAGS=-O2 -std=c++14

#CXXFLAGS+=-O0 -g

HEADS=ntt.h mod.h base.h
HEADS+=ref.h
#HEADS+=pmp.h
#HEADS+=dpmp.h
#HEADS+=pmp2.h
#HEADS+=pmp5.h
#HEADS+=cooley.h
#HEADS+=stockham_dit.h
#HEADS+=stockham_dif.h
#HEADS+=stockham6.h

.PHONY: all
.PHONY: run

all: run

run: ntt_bench
	./ntt_bench

ntt_bench: main.cc $(HEADS)
	$(CXX) $(CXXFLAGS) $< -o $@

CC=clang
CXX=clang++
CCFLAGS=-O2
CXXFLAGS=$(CCFLAGS) -std=c++17


all: itchyny peria


itchyny: itchyny.c
	$(CC) -o $@ $^ $(CCFLAGS) -lgmp -lm

peria: chudnovsky.cc computer.cc peria.cc ramanujan.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) -lgmpxx -lgmp

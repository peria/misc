CC=clang
CXX=clang++
CCFLAGS=-O2
CXXFLAGS=$(CCFLAGS) -std=c++17


all: itchyny peria


itchyny: itchyny.c
	$(CC) -o $@ $^ $(CCFLAGS) -lgmp -lm

peria: peria.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) -lgmpxx -lgmp

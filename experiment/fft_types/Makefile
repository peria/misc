CXX=clang++
CXXFLAGS=-O2 -std=c++14

# CXXFLAGS+=-O0 -g

SRCS=main.cc fmt.cc
HEADS=fmt.h complex.h
HEADS+=dit.h
OBJS=$(SRCS:.cc=.o)

all: run

run: benchmark
	./$<

benchmark: $(OBJS) $(HEADS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

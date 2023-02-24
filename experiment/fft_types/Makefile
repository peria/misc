CXX=clang++
CXXFLAGS=-O2 -std=c++14

CXXFLAGS+=-O0 -g

SRCS=main.cc integer.cc fmt.cc
HEADS=fmt.h complex.h integer.h
# HEADS+=pmp.h
HEADS+=dpmp.h
HEADS+=ct_dit.h
# HEADS+=pmp2.h
# HEADS+=pmp5.h
# HEADS+=cooley.h
# HEADS+=stockham_dit.h
# HEADS+=stockham_dif.h
# HEADS+=stockham6.h
OBJS=$(SRCS:.cc=.o)

all: run

run: benchmark
	./$<

benchmark: $(OBJS) $(HEADS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

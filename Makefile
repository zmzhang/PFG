CC          = gcc
CXX         = g++
SHELL       = sh
AR = ar cr
RANLIB = ranlib

CFLAGS = -Wall -O3
LDFLAGS = 

ifneq ($(NO_OPENMP), 1)
CFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif

SRCS = compound.cpp elements.cpp formula.cpp getopt_pp.cpp
OBJS = $(SRCS:.cpp=.o)

PFG: $(OBJS)
	$(CXX) -o PFG $(OBJS) $(LDFLAGS)

.PHONY: clean

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

all: PFG
	
clean:
	rm -f PFG
	rm -f PFG.exe
	rm -f *.o
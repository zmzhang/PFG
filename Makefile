CC          = gcc
CXX         = g++
SHELL       = sh
AR = ar cr
RANLIB = ranlib

PFG: *.cpp
	$(CXX) -o PFG meta.h compound.cpp elements.cpp formula.cpp getopt_pp.cpp -O3 -fopenmp 
	
clean:
	rm -f PFG.exe
	rm -f *.o
#include *.mk

all : fzip

fzip : fzip.cpp
	g++ -o fzip fzip.cpp

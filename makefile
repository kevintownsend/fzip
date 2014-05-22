include *.mk

all : fzip

fzip : src/fzip.cpp
	g++ -o fzip src/fzip.cpp

all : fzip

fzip : source/fzip.cpp source/arithmeticFloatingPointCompression.o
	g++ -std=c++0x -Iinclude -o fzip source/fzip.cpp source/*.o

source/arithmeticFloatingPointCompression.o : source/arithmeticFloatingPointCompression.cpp
	g++ -Iinclude -c -std=c++0x source/arithmeticFloatingPointCompression.cpp -o source/arithmeticFloatingPointCompression.o

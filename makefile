all : fzip

fzip : source/fzip.cpp source/arithmeticFloatingPointCompression.o source/common.o source/sortByRepeats.o
	g++ -std=gnu++0x -L${HOME}/lib -ltardis -Iinclude -o fzip source/fzip.cpp source/*.o

source/arithmeticFloatingPointCompression.o : source/arithmeticFloatingPointCompression.cpp
	g++ -Iinclude -L${HOME}/lib -ltardis -I${HOME}/include -Iinclude -c -std=c++0x source/arithmeticFloatingPointCompression.cpp -o source/arithmeticFloatingPointCompression.o

source/common.o : source/common.cpp
	g++ -Iinclude -c -std=gnu++0x source/common.cpp -o source/common.o

source/sortByRepeats.o : source/sortByRepeats.cpp
	g++ -Iinclude -c -std=gnu++0x source/sortByRepeats.cpp -o source/sortByRepeats.o

install :
	cp fzip ${HOME}/bin/.

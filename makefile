all: fzip

run :
	fzip

fzip: source/fzip.cpp source/common.o source/sortByRepeats.o include/fzip.hpp
	g++ -std=c++11 -Iinclude -o fzip source/fzip.cpp

source/arithmeticFloatingPointCompression.o: source/arithmeticFloatingPointCompression.cpp
	g++ -Iinclude -L${HOME}/lib -ltardis -I${HOME}/include -Iinclude -c -std=c++0x source/arithmeticFloatingPointCompression.cpp -o source/arithmeticFloatingPointCompression.o

source/common.o: source/common.cpp
	g++ -Iinclude -c -std=gnu++0x source/common.cpp -o source/common.o

source/sortByRepeats.o: source/sortByRepeats.cpp
	g++ -Iinclude -c -std=gnu++0x source/sortByRepeats.cpp -o source/sortByRepeats.o

source/fzipCore.o: source/fzipCore.cpp
	g++ -Iinclude -c -std=c++11 source/fzipCore.cpp -o source/fzipCore.o

install:
	cp fzip ${HOME}/bin/.

clean:
	rm -rf fzip source/*.o

vim:
	vim -p makefile ./include/fzip.hpp ./source/fzip.cpp ../spMatrixHelp/rcr.hpp ./source/fzipCore.cpp

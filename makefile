all: fzip floatDiff

run : fzip
	fzip -c < ./benchmark/example.bin > example.fz
	fzip -d example.fz exampleAfter.bin

run2 : fzip
	fzip -c ./benchmark/obs_info.trace.bin obs_info.fz
	fzip -d obs_info.fz obs_info.after

run3: fzip
	fzip -c ./benchmark/test.trace.bin test.fz
	fzip -d test.fz test.after


fzip: source/fzip.cpp source/common.o source/sortByRepeats.o include/fzip.hpp
	g++ -g -std=c++11 -Iinclude -o fzip source/fzip.cpp

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


floatDiff: floatDiff.cpp
	g++ -o floatDiff floatDiff.cpp

vim:
	vim -p makefile ./include/fzip.hpp ./source/fzip.cpp floatDiff.cpp ../spMatrixHelp/rcr.hpp ./source/fzipCore.cpp

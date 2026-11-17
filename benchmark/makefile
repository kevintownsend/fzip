CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -I../include

.PHONY: all clean

all:
	$(MAKE) -C .. bench

clean:
	rm -f *.fz *.bin bench

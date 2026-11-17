CXX ?= g++
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCDIR ?= $(PREFIX)/include
LIBDIR ?= $(PREFIX)/lib

CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -Wpedantic -Werror -Iinclude
ifeq ($(SANITIZE),1)
OPTFLAGS := -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS_EXTRA := -fsanitize=address,undefined
else ifeq ($(DEBUG),1)
OPTFLAGS := -g -O1
LDFLAGS_EXTRA :=
else
OPTFLAGS := -O3
LDFLAGS_EXTRA :=
endif

CXXFLAGS := $(CXXFLAGS_COMMON) $(OPTFLAGS) $(EXTRA_CXXFLAGS)

BUILD_DIR := build
LIB_OBJS := $(BUILD_DIR)/fzip.o
CLI_OBJS := $(BUILD_DIR)/main.o
TEST_SRCS := $(filter-out tests/fuzz_%,$(wildcard tests/*.cpp))
BENCH_SRCS := $(wildcard benchmark/*.cpp)

.PHONY: all lib test bench install clean

all: fzip

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/fzip.o: src/fzip.cpp include/fzip/fzip.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/main.o: src/main.cpp include/fzip/fzip.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/libfzip.a: $(LIB_OBJS)
	$(AR) rcs $@ $^

lib: $(BUILD_DIR)/libfzip.a

fzip: $(CLI_OBJS) $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS_EXTRA) $^ -o $@

test: fzip
ifneq ($(TEST_SRCS),)
	$(CXX) $(CXXFLAGS) $(LDFLAGS_EXTRA) $(TEST_SRCS) $(LIB_OBJS) -o $(BUILD_DIR)/fzip_test
	./$(BUILD_DIR)/fzip_test
else
	./fzip --help >/dev/null
	@printf '\0\0\0\0\0\0\0\0\0\0\0\0\0\0\360\77' | ./fzip -c | ./fzip -d >/dev/null
	@echo "CLI smoke test passed (add tests/*.cpp for full suite)"
endif

bench: $(LIB_OBJS)
ifneq ($(BENCH_SRCS),)
	$(CXX) $(CXXFLAGS) $(LDFLAGS_EXTRA) $(BENCH_SRCS) $(LIB_OBJS) -o $(BUILD_DIR)/bench
	./$(BUILD_DIR)/bench
else
	@echo "No benchmark sources found under benchmark/"
endif

install: fzip $(BUILD_DIR)/libfzip.a
	install -d $(DESTDIR)$(BINDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(INCDIR)/fzip
	install -m 755 fzip $(DESTDIR)$(BINDIR)/fzip
	install -m 644 $(BUILD_DIR)/libfzip.a $(DESTDIR)$(LIBDIR)/libfzip.a
	install -m 644 include/fzip/fzip.hpp $(DESTDIR)$(INCDIR)/fzip/fzip.hpp

clean:
	rm -rf $(BUILD_DIR) fzip
	$(MAKE) -C benchmark clean

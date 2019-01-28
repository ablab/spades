# Makefile for suftest and test

# options
CC						= llvm-gcc
#CXX						= g++
#OUTPUT_OPTION	= -o $@
CFLAGS                          = -ffast-math -O9 -funroll-loops -DNDEBUG
#CFLAGS				= -O3 -fomit-frame-pointer -funroll-loops
#CXXFLAGS			= -O3 -fomit-frame-pointer
CPPFLAGS			= -Wall -DNDEBUG
#CPPFLAGS			= -Wall
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=

# targets
.PHONY: all
all: suftest
suftest: sais.o suftest.o
test:
	$(CC) -O -g -Wall test.c sais.c -o test
	./test
	$(RM) test test.exe

distclean: clean
clean:
	$(RM) suftest suftest.exe test test.exe sais.o suftest.o

# dependencies
sais.o suftest.o: sais.h Makefile

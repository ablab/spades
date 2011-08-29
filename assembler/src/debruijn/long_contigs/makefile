OBJ_DIR=../../../build/release/debruijn/CMakeFiles/debruijn.dir
MYO_DIR=../../../build/release/debruijn/long_contigs
APP_DIR=../../../build/release/debruijn/long_contigs

CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=gnu++0x $(PG) $(OMP) -I./../../include/io  -I./../../../ext/include  -I./../../../ext/include/statgen -I./../../include/read -I./../../sequence -I./../../include/omni -I./../../include/libs -I./.. -I./../../include

LDFLAGS=-L../../../build/release/io -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -linput 

CC=g++

default: long_contigs

long_contigs: $(MYO_DIR)/main.o
	$(CC) $(LDFLAGS) $(MYO_DIR)/main.o  -o $(APP_DIR)/long_contigs

$(MYO_DIR)/main.o: main.cpp lc_config_struct.hpp lc_common.hpp seeds.hpp paths.hpp quality.hpp visualize.hpp path_utils.hpp loop.hpp  lc_io.hpp
	mkdir -p $(MYO_DIR)
	$(CC) $(CFLAGS) main.cpp -o $(MYO_DIR)/main.o

clean: 
	rm $(MYO_DIR)/*.o -rf

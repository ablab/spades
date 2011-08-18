OBJ_DIR=../../../build/release/debruijn/CMakeFiles/debruijn.dir
MYO_DIR=../../../build/release/debruijn/long_contigs
APP_DIR=../../../build/release/debruijn/long_contigs

CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=gnu++0x $(PG) $(OMP) -I./../../common -I./../../common/read -I./../../common/sequence -I./../../omni -I./../../libs -I./.. -I./../..

LDFLAGS=-O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1

CC=g++

default: long_contigs

long_contigs: main.o 
	$(CC) $(LDFLAGS) $(MYO_DIR)/main.o  -o $(APP_DIR)/long_contigs

main.o: main.cpp lc_config_struct.hpp lc_common.hpp seeds.hpp paths.hpp quality.hpp visualize.hpp ../new_debruijn.hpp ../../omni/visualization_utils.hpp
	$(CC) $(CFLAGS) main.cpp -o $(MYO_DIR)/main.o

clean: 
	rm $(MYO_DIR)/*.o -rf



DIR=../../build/tools/cuckoo_test

program: cuckoo_test.o sequence.o
	g++ $(DIR)/cuckoo_test.o $(DIR)/sequence.o -O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1 -o $(DIR)/cuckoo_test

cuckoo_test.o:
	mkdir -p $(DIR)
	rm -rf $(DIR)/*
	g++ -c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=c++0x -I./../../src/common -I./../../src -I./.. -I./../../src/common/read -I./../../src/common/sequence cuckoo_test.cpp -o ../../build/tools/cuckoo_test/cuckoo_test.o

sequence.o:
	g++ -c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=c++0x -I./../../src/common -I./../../src  -I./../../src/common/read -I./../../src/common/sequence ../../src/common/sequence/sequence.cpp -o ../../build/tools/cuckoo_test/sequence.o

clean:
	rm -rf ../../build/tools/cuckoo_test/*

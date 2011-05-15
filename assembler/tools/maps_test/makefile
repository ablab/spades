DIR=../../build/tools/maps_test

program: filter.o sequence.o
	g++ $(DIR)/filter.o $(DIR)/sequence.o -O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1 -o $(DIR)/filter

filter.o: 
	mkdir -p $(DIR)
	rm -rf $(DIR)/*
	g++ -c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=c++0x   -I./../../src/common -I./../../src -I./..  filter.cpp -o $(DIR)/filter.o

sequence.o:
	g++ -c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=c++0x   -I./../../src/common -I./../../src   ../../src/common/sequence.cpp -o $(DIR)/sequence.o

clean:
	rm -rf $(DIR)/*

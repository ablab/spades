#include <iostream>
#include <cstdlib>
#include "seq_filter.hpp"

#define K 31

int main(int argc, char** argv) {
	if (argc == 3) {
		std::string filename(argv[1]);
		size_t L = atoi(argv[2]);
		seq_filter<K>::filter(filename, L);
	} else {
		std::cout << "Usage: ./filter <filename> <L>\n";
		std::cout << "Selects k-mer with amount > L\n";
	}
  return 0;
}

#include <iostream>
#include <cstdlib>
#include "seq_filter.hpp"

int main(int argc, char** argv) {
	if (argc == 4) {
		std::string filename(argv[1]);
		const size_t size = 5; //atoi(argv[2]);
		size_t L = atoi(argv[3]);
		seq_filter<size>::filter(filename, L);
	} else {
		std::cout << "Usage: ./filter <filename> <size> <L>\n";
	}
  return 0;
}

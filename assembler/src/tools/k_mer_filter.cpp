#include <iostream>
#include "seq_filter.hpp"

const std::string filename = "../../src/tools/s_6.first10000_1.fastq.gz";
//const std::string filename = "/dev/tty";
//const std::string output = "/dev/stdout";
const size_t size_ = 5;

int main() {
	seq_filter<size_>::filter(filename, 5);
  return 0;
}

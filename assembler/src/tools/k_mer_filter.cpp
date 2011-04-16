//#include "seq.hpp"
#include "strobe_read.hpp"
//#include "nucl.hpp"
//#include "sequence.hpp"
//#include "read.hpp"
#include "ireadstream.hpp"
#include <iostream>
#include "cuckoo.hpp"

std::string filename = "../../src/tools/s_6.first10000_1.fastq.gz";

int main() {
  //ireadstream ifs(filename);
  //std::vector<Read> reads;
  //Read r;
  //strobe_read<10> kmer;
  //while (!ifs.eof()) {
  //  ifs >> r;
  //  reads.push_back(r);
  //}
  //ifs.close();
  vector<Read>* reads = ireadstream::readAll(filename, 10000);

  for (size_t i = 0; i < reads->size(); i++) {
    std::string s = (*reads)[i].getSequenceString();
    strobe_read<5, 1> kmer(&s);
    std::cout << s << " " << s.size() << " " << kmer[0] << std::endl;
  }
  std::cout << reads->size() << std::endl;

  return 0;
}

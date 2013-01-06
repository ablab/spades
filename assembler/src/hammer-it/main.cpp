#include "logger/log_writers.hpp"

#include "read/read.hpp"
#include "read/ireadstream.hpp"

#include "HSeq.hpp"
#include "kmer_data.hpp"
#include "valid_hkmer_generator.hpp"

void create_console_logger() {
	using namespace logging;

	logger *lg = create_logger("");
	lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}
 
int main(void) {
  srand(42);
  srandom(42);

  create_console_logger();

  KMerData km;
  KMerDataCounter(1).FillKMerData(km);
  
  return 0;
}

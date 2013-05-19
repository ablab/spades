#include "running_modes.hpp"

#include "adapter_index.hpp"
#include "job_wrappers.hpp"
#include "output.hpp"

#include "io/read_processor.hpp"
#include "logger/log_writers.hpp"


void exactAndAlign(std::ostream& output, std::ostream& bed, ireadstream * input,
                   const cclean::AdapterIndex &index) {
  ExactAndAlignJobWrapper filler(output, bed, index);
  hammer::ReadProcessor rp(cfg::get().nthreads);
  rp.Run(*input, filler);
  VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
  INFO("Reads processed: " << rp.processed() << ". Reads aligned: " << filler.aligned());
}

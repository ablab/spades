#include "running_modes.hpp"

#include <string>

#include "adapter_index.hpp"
#include "job_wrappers.hpp"
#include "output.hpp"
#include "io/read_processor.hpp"
#include "logger/log_writers.hpp"
#include "additional.cpp"
#include "brute_force_clean.hpp"

using additional::WorkModeType;
using additional::NONE;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;
using std::string;

void ExactAndAlign(std::ostream &output, std::ostream &bed, ireadstream *input,
                   const string &db, const cclean::AdapterIndex &index,
                   const WorkModeType &mode) {
  if (mode == SIMPLE) {
    SimpleClean filler(output, bed, db, index);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
  }
  else {  // Now only two modes, thats because else
    BruteForceClean filler(output, bed, db, index.GetSeqs(), mode);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
  }
}

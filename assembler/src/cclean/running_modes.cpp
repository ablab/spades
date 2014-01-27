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
using additional::SIGNLE_END;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_WITH_Q;
using std::string;

void ExactAndAlign(std::ostream &aligned_output,
                   std::ostream &bed, ireadstream *input, std::ostream &output,
                   const std::string &db, const cclean::AdapterIndex &index,
                   const additional::WorkModeType &mode) {
  if (mode == SIGNLE_END) {
    SimpleClean filler(aligned_output, output, bed, db, index);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
    INFO("Reads aligned: " << filler.aligned());
  }
  else {  // Now only two modes, thats because else
    BruteForceClean filler(aligned_output, output, bed, db, index.GetSeqs(), mode);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
  }
}

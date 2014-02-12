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
using additional::SINGLE_END;
using additional::SINGLE_END_Q;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_WITH_Q;
using std::string;

void ExactAndAlign(std::ostream &aligned_output,
                   std::ostream &bed, ireadstream *input, std::ostream &output,
                   const std::string &db, const cclean::AdapterIndex &index,
                   const additional::WorkModeType &mode,
                   const std::unordered_map<string, string> &options) {
  if (mode == SINGLE_END || mode == SINGLE_END_Q) {
    SimpleClean filler(aligned_output, output, bed, db, index, mode, options);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
    INFO("Reads aligned: " << filler.aligned());
  }
  else {
    BruteForceClean filler(aligned_output, output, bed, db, index.GetSeqs(),
                           mode, options);
    hammer::ReadProcessor rp(cfg::get().nthreads);
    rp.Run(*input, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
    INFO("Reads processed: " << rp.processed() << ". Reads aligned: "
         << filler.aligned());
  }
}

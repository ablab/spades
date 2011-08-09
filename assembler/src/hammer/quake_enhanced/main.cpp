#include <iostream>
#include "common/logging.hpp"
#include "hammer/quake_enhanced/options.hpp"
#include "hammer/quake_enhanced/quake.hpp"
DECL_PROJECT_LOGGER("q")

namespace quake_enhanced {

DECL_LOGGER("main")

}

int main(int argc, char **argv) {
  quake_enhanced::Options opts(argc, argv);
  if (!opts.valid) {
    std::cout << opts.help_message;
    return 1;
  }
  quake_enhanced::Quake quake;
  quake.Count(opts.read_file, opts.kmer_count_file,
              opts.hash_file_prefix, opts.hash_file_number,
              opts.quality_offset, opts.quality_threshold);
  quake.PrepareHists(opts.hist_file, opts.trusted_hist_file, 
                     opts.bad_hist_file, opts.top_threshold,
                     opts.average_min);
  return 0;
}

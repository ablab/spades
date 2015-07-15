//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include "logging.hpp"
#include "options.hpp"
#include "quake.hpp"


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
  quake.PrepareLimits(opts.bad_threshold, opts.limits_file);
  quake.FilterTrusted(opts.kmer_count_file, opts.trusted_hist_file, opts.bad_kmer_file);
  return 0;
}

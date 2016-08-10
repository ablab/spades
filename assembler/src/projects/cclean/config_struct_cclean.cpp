//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "config_struct_cclean.hpp"
#include "pipeline/config_common.hpp"
#include "utils/openmp_wrapper.h"

void load(cclean_config& cfg, const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  load(cfg, pt);
}

void load(cclean_config& cfg, boost::property_tree::ptree const& pt) {
  using config_common::load;
  load(cfg.use_quality, pt, "use_quality");
  load(cfg.use_bruteforce, pt, "use_bruteforce");
  load(cfg.debug_information, pt, "debug_information");

  load(cfg.score_treshold, pt, "score_treshold");
  load(cfg.mismatch_threshold, pt, "mismatch_threshold");
  load(cfg.minimum_lenght, pt, "minimum_lenght");
  load(cfg.nthreads, pt, "nthreads");
  load(cfg.aligned_part_fraction, pt, "aligned_part_fraction");
  load(cfg.buffer_size, pt, "buffer_size");

  load(cfg.dataset_file_name, pt, "dataset");
  load(cfg.database, pt, "database");
  load(cfg.input_working_dir, pt, "input_working_dir");
  load(cfg.output_working_dir, pt, "output_working_dir");

  std::string file_name = cfg.dataset_file_name;
  cfg.dataset.load(file_name);

  // Fix number of threads according to OMP capabilities.
  cfg.nthreads = std::min(cfg.nthreads, (unsigned)omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads(cfg.nthreads);
}

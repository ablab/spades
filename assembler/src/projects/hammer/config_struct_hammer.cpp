//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * config_struct_hammer.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: snikolenko
 */

#include "config_struct_hammer.hpp"
#include "pipeline/config_common.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <boost/property_tree/ptree.hpp>
#include <string>

void load(hammer_config& cfg, const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  load(cfg, pt);
}

void load(hammer_config& cfg, boost::property_tree::ptree const& pt) {
  using config_common::load;
  load(cfg.general_do_everything_after_first_iteration, pt, "general_do_everything_after_first_iteration");
  load(cfg.general_hard_memory_limit, pt, "general_hard_memory_limit");
  load(cfg.general_max_nthreads, pt, "general_max_nthreads");
  load(cfg.general_tau, pt, "general_tau");
  load(cfg.general_max_iterations, pt, "general_max_iterations");
  load(cfg.general_debug, pt, "general_debug");

  load(cfg.count_do, pt, "count_do");
  load(cfg.count_numfiles, pt, "count_numfiles");
  load(cfg.count_merge_nthreads, pt, "count_merge_nthreads");
  load(cfg.count_split_buffer, pt, "count_split_buffer");
  load(cfg.count_filter_singletons, pt, "count_filter_singletons");
  
  load(cfg.hamming_do, pt, "hamming_do");
  load(cfg.hamming_blocksize_quadratic_threshold, pt, "hamming_blocksize_quadratic_threshold");

  load(cfg.bayes_do, pt, "bayes_do");
  load(cfg.bayes_nthreads, pt, "bayes_nthreads");
  load(cfg.bayes_singleton_threshold, pt, "bayes_singleton_threshold");
  load(cfg.bayes_nonsingleton_threshold, pt, "bayes_nonsingleton_threshold");
  load(cfg.bayes_discard_only_singletons, pt, "bayes_discard_only_singletons");
  load(cfg.bayes_debug_output, pt, "bayes_debug_output");
  load(cfg.bayes_use_hamming_dist, pt, "bayes_use_hamming_dist");
  load(cfg.bayes_hammer_mode, pt, "bayes_hammer_mode");
  load(cfg.bayes_write_bad_kmers, pt, "bayes_write_bad_kmers");
  load(cfg.bayes_write_solid_kmers, pt, "bayes_write_solid_kmers");
  load(cfg.bayes_initial_refine, pt, "bayes_initial_refine");

  load(cfg.expand_do, pt, "expand_do");
  load(cfg.expand_max_iterations, pt, "expand_max_iterations");
  load(cfg.expand_nthreads, pt, "expand_nthreads");
  load(cfg.expand_write_each_iteration, pt, "expand_write_each_iteration");
  load(cfg.expand_write_kmers_result, pt, "expand_write_kmers_result");

  load(cfg.correct_do, pt, "correct_do");
  load(cfg.correct_nthreads, pt, "correct_nthreads");
  load(cfg.correct_threshold, pt, "correct_threshold");
  load(cfg.correct_use_threshold, pt, "correct_use_threshold");
  load(cfg.correct_readbuffer, pt, "correct_readbuffer");
  load(cfg.correct_discard_bad, pt, "correct_discard_bad");
  load(cfg.correct_stats, pt, "correct_stats");

  std::string fname;
  load(fname, pt, "dataset");
  cfg.dataset.load(fname);

  load(cfg.input_working_dir, pt, "input_working_dir");
  load(cfg.input_trim_quality, pt, "input_trim_quality");
  cfg.input_qvoffset_opt = pt.get_optional<int>("input_qvoffset");
  load(cfg.output_dir, pt, "output_dir");

  cfg.general_max_nthreads = spades_set_omp_threads(cfg.general_max_nthreads);
}

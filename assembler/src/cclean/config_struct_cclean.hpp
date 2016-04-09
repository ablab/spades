//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef CONFIG_STRUCT_CCLEAN_HPP
#define CONFIG_STRUCT_CCLEAN_HPP

#include "pipeline/config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include "pipeline/library.hpp"

struct cclean_config {

  bool use_quality;
  bool use_bruteforce;
  bool debug_information;

  unsigned score_treshold;
  unsigned mismatch_threshold;
  unsigned minimum_lenght;
  unsigned nthreads;
  unsigned buffer_size;
  double aligned_part_fraction;

  std::string dataset_file_name;
  std::string database;
  std::string input_working_dir;
  std::string output_working_dir;

  io::DataSet<> dataset;
};

// main config load function
void load(cclean_config& cfg, const std::string &filename);
void load(cclean_config& cfg, boost::property_tree::ptree const& pt);

typedef config_common::config<cclean_config> cfg;

#endif

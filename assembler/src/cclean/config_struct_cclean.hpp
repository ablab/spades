//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef CONFIG_STRUCT_CCLEAN_HPP_
#define CONFIG_STRUCT_CCLEAN_HPP_

#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct cclean_config {
  unsigned mismatch_threshold;
  double aligned_part_fraction;
  std::string output_file;
  std::string output_bed;
  unsigned nthreads;

  std::string input_working_dir;
};

// main config load function
void load(cclean_config& cfg, const std::string &filename);
void load(cclean_config& cfg, boost::property_tree::ptree const& pt);

typedef config_common::config<cclean_config> cfg;

#endif

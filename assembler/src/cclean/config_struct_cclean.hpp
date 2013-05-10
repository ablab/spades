//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * config_struct_cclean.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: snikolenko
 */

#ifndef CONFIG_STRUCT_CCLEAN_HPP_
#define CONFIG_STRUCT_CCLEAN_HPP_

#include "config_common.hpp"
#include "boost/optional.hpp"

// struct for debruijn project's configuration file
struct cclean_config
{
  int mismatch_threshold;
  double aligned_part_fraction;
  std::string output_file;
  std::string output_bed;
  int nthreads;
};


// main config load function
void load(cclean_config& cfg, boost::property_tree::ptree const& pt);

typedef config_common::config<cclean_config> cclean_cfg;

#endif

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_common.hpp"
#include "pipeline/config_singl.hpp"

namespace cap {

struct cap_config {
  std::string cache_root;
  std::string desc_file_name;
  std::string default_log_filename;
  std::string default_log_file_mode;
};

inline void load(cap_config &cfg, boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;

  load(cfg.cache_root, pt, "cache_root");
  load(cfg.desc_file_name, pt, "desc_file_name");
  load(cfg.default_log_filename, pt, "default_log_filename");
  load(cfg.default_log_file_mode, pt, "default_log_file_mode");
}

void load(cap_config& cfg, const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  load(cfg, pt, true);
}

}

typedef config_common::config<cap::cap_config> cap_cfg;

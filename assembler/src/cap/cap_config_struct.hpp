#pragma once

#include "config_common.hpp"
#include <vector>

namespace cap {

struct cap_config {
  std::string cache_root;
  std::string desc_file_name;
};

inline void load(cap_config &cfg, boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

  load(cfg.cache_root, pt, "cache_root");
  load(cfg.desc_file_name, pt, "desc_file_name");
}

}

typedef config_common::config<cap::cap_config> cap_cfg;

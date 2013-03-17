#ifndef __HAMMER_IT_CONFIG_HPP__
#define __HAMMER_IT_CONFIG_HPP__

#include "config_singl.hpp"

struct hammer_config {
};

void load(hammer_config& cfg, const std::string &filename);

typedef config_common::config<hammer_config> cfg;

#endif // __HAMMER_IT_CONFIG_HPP__

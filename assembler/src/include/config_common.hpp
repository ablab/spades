//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * config_common.hpp
 *
 *  Created on: Aug 13, 2011
 *      Author: Alexey.Gurevich
 */

#pragma once

#include "simple_tools.hpp"
#include "path_helper.hpp"
#include "verify.hpp"

// todo: undo dirty fix
#include <boost/format.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>

namespace config_common {
// for enable_if/disable_if
namespace details {
template<class T, class S>
struct is_equal_type {
  static const bool value = false;
};

template<class T>
struct is_equal_type<T, T> {
  static const bool value = true;
};
}

template<class T>
typename boost::enable_if_c<
  details::is_equal_type<T, std::string>::value ||
  boost::is_arithmetic<T>::value>::type load(T& value,
                                             boost::property_tree::ptree const& pt, std::string const& key,
                                             bool complete) {
  if (complete || pt.find(key) != pt.not_found())
    value = pt.get<T>(key);
}

template<class T>
typename boost::disable_if_c<
  details::is_equal_type<T, std::string>::value ||
  boost::is_arithmetic<T>::value>::type load(T& value,
                                             boost::property_tree::ptree const& pt, std::string const& key,
                                             bool complete) {
  if (complete || pt.find(key) != pt.not_found())
    load(value, pt.get_child(key), complete);
}

template<class T>
void load_items(std::vector<T>& vec, boost::property_tree::ptree const& pt,
                std::string const& key, bool complete) {
  std::string vector_key = key + std::string(".count");
  if (complete || pt.find(vector_key) != pt.not_found()) {
    size_t count = pt.get<size_t>(vector_key);

    for (size_t i = 0; i != count; ++i) {
      T t;
      load(t, pt.get_child(str(boost::format("%s.item_%d") % key % i)),
           complete);
      vec.push_back(t);
    }
  }
}

void inline split(std::vector<std::string>& vec, std::string const& space_separated_list) {
  std::istringstream iss(space_separated_list);
  while (iss) {
    std::string value;
    iss >> value;
    if (value.length()) {
      vec.push_back(value);
    }
  }
}

void inline load_split(std::vector<std::string>& vec, boost::property_tree::ptree const& pt, std::string const& key) {
  boost::optional<std::string> values = pt.get_optional<std::string>(key);
  if (values) {
    split(vec, *values);
  }
}

template<class T>
void inline load(std::vector<T>& vec, boost::property_tree::ptree const& pt, std::string const& key, bool /*complete*/) {
	boost::optional<T> value = pt.get_optional<T>(key);
	if (value) {
		vec.push_back(*value);
		return;
	}
	for (size_t i = 1;; i++) {
		value = pt.get_optional<std::string>(key + "#" + ToString(i));
		if (value) {
			vec.push_back(*value);
			continue;
		}
		value = pt.get_optional<std::string>(key + "." + ToString(i));
		if (value) {
			vec.push_back(*value);
			continue;
		}
		if (i > 0) {
			return;
		}
	}
}

template<class T>
void load(T& value, boost::property_tree::ptree const& pt, std::string const& key) {
  load(value, pt, key, true);
}

template<class T>
void load(T& value, boost::property_tree::ptree const& pt, const char* key) {
  load(value, pt, std::string(key), true);
}

template<class T>
void load(T& value, boost::property_tree::ptree const& pt) {
  load(value, pt, true);
}
}

template<class T>
inline void load_param(const std::string& filename, const std::string& key,
                       boost::optional<T>& value) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);
  value = pt.get_optional<T>(key);
}

template<class T>
inline void write_param(const std::string& filename, const std::string& key,
                        const boost::optional<T>& value) {
  if (value) {
    std::ofstream params_stream(filename.c_str(), std::ios_base::app);
    params_stream << key << "\t" << value << std::endl;
  }
}

template<class T>
inline void write_param(const std::string& filename, const std::string& key,
                        const T &value) {
  std::ofstream params_stream(filename.c_str(), std::ios_base::app);
  params_stream << key << "\t" << value << std::endl;
}

template<class K, class V>
inline void load_param_map(const std::string& filename, const std::string& key,
                           std::map<K, V>& value) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);
  boost::optional<std::string> as_str = pt.get_optional<std::string>(key);
  if (as_str) {
    std::vector<std::string> key_value_pairs;
    boost::split(key_value_pairs, *as_str, boost::is_any_of(";"));
    for (auto it = key_value_pairs.begin(); it != key_value_pairs.end();
         ++it) {
      std::vector<std::string> key_value;
      boost::split(key_value, *it, boost::is_any_of(" "));
      VERIFY(key_value.size() == 2);
      value[boost::lexical_cast<K>(key_value[0])] =
          boost::lexical_cast<K>(key_value[1]);
    }
  }
}

template<class K, class V>
inline void write_param_map(const std::string& filename, const std::string& key,
                            const std::map<K, V>& value) {
  if (value.size() > 0) {
    std::ofstream params_stream(filename.c_str(), std::ios_base::app);
    params_stream << key << "\t\"";
    std::string delim = "";
    for (auto it = value.begin(); it != value.end(); ++it) {
      params_stream << delim << it->first << " " << it->second;
      delim = ";";
    }
    params_stream << "\"" << std::endl;
  }
}

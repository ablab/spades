//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/stl_utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/verify.hpp"

// todo: undo dirty fix

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

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
typename boost::enable_if_c<details::is_equal_type<T, std::string>::value ||
                            boost::is_arithmetic<T>::value>::type
load(T &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found())
        value = pt.get<T>(key);
}

template<class T>
typename boost::disable_if_c<details::is_equal_type<T,
                                                    std::string>::value ||
                             boost::is_arithmetic<T>::value>::type
load(T &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found())
        load(value, pt.get_child(key), complete);
}

template<class T>
void load_items(std::vector <T> &vec, boost::property_tree::ptree const &pt,
                std::string const &key, bool complete) {
    std::string vector_key = key + std::string(".count");
    if (complete || pt.find(vector_key) != pt.not_found()) {
        size_t count = pt.get<size_t>(vector_key);

        for (size_t i = 0; i != count; ++i) {
            T t;
            load(t, pt.get_child(fmt::format("{:s}.item_{:d}", key, i)),
                 complete);
            vec.push_back(t);
        }
    }
}

template<class T>
void load(std::vector <T> &vec, boost::property_tree::ptree const &pt, std::string const &key,
          bool /*complete*/) {
    boost::optional<T> value = pt.get_optional<T>(key);
    if (value) {
        vec.push_back(*value);
        return;
    }
    for (size_t i = 1; ; i++) {
        value = pt.get_optional<std::string>(key + "#" + std::to_string(i));
        if (value) {
            vec.push_back(*value);
            continue;
        }
        value = pt.get_optional<std::string>(key + "." + std::to_string(i));
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
void load(T &value, boost::property_tree::ptree const &pt, std::string const &key) {
    load(value, pt, key, true);
}

template<class T>
void load(T &value, boost::property_tree::ptree const &pt, const char *key) {
    load(value, pt, std::string(key), true);
}

template<class T>
void load(T &value, boost::property_tree::ptree const &pt) {
    load(value, pt, true);
}

template<class T>
void load_param(const std::string &filename, const std::string &key,
                boost::optional<T> &value) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    value = pt.get_optional<T>(key);
}

template<class T>
void write_param(const std::string &filename, const std::string &key,
                 const boost::optional<T> &value) {
    if (value) {
        std::ofstream params_stream(filename.c_str(), std::ios_base::app);
        params_stream << key << "\t" << value << std::endl;
    }
}

template<class T>
void write_param(const std::string &filename, const std::string &key,
                 const T &value) {
    std::ofstream params_stream(filename.c_str(), std::ios_base::app);
    params_stream << key << "\t" << value << std::endl;
}

}

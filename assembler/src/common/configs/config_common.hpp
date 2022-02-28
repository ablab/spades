//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/stl_utils.hpp"
#include "utils/verify.hpp"

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <cppformat/format.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

namespace config_common {

template<class T>
typename boost::enable_if_c<std::is_convertible<T, std::string>::value ||
                            boost::is_arithmetic<T>::value>::type
load(T &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found())
        value = pt.get<T>(key);
}

template<class T>
typename boost::enable_if_c<std::is_convertible<T, std::string>::value ||
                            boost::is_arithmetic<T>::value>::type
load(std::optional<T> &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found())
        value = pt.get<T>(key);
}

template<class T>
typename boost::disable_if_c<std::is_convertible<T, std::string>::value ||
                             boost::is_arithmetic<T>::value>::type
load(T &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found())
        load(value, pt.get_child(key), complete);
}

template<class T>
typename boost::disable_if_c<std::is_convertible<T, std::string>::value ||
                             boost::is_arithmetic<T>::value>::type
load(std::optional<T> &value,
     boost::property_tree::ptree const &pt, std::string const &key,
     bool complete) {
    if (complete || pt.find(key) != pt.not_found()) {
        T new_value;
        load(new_value, pt.get_child(key), complete);
        value = new_value;
    }
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
    std::optional<T> value = pt.get_optional<T>(key);
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
void load_param(const std::filesystem::path &filename, const std::string &key,
                std::optional<T> &value) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    value = pt.get_optional<T>(key);
}

template<class T>
void write_param(const std::filesystem::path &filename, const std::string &key,
                 const std::optional<T> &value) {
    if (value) {
        std::ofstream params_stream(filename, std::ios_base::app);
        params_stream << key << "\t" << value << std::endl;
    }
}

template<class T>
void write_param(const std::filesystem::path &filename, const std::string &key,
                 const T &value) {
    std::ofstream params_stream(filename, std::ios_base::app);
    params_stream << key << "\t" << value << std::endl;
}

}

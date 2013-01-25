//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

// todo: undo dirty fix
#include <boost/format.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "simple_tools.hpp"
#include "path_helper.hpp"
#include "verify.hpp"

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
		details::is_equal_type<T, std::string>::value
				|| boost::is_arithmetic<T>::value>::type load(T& value,
		boost::property_tree::ptree const& pt, string const& key,
		bool complete) {
	if (complete || pt.find(key) != pt.not_found())
		value = pt.get<T>(key);
}

template<class T>
typename boost::disable_if_c<
		details::is_equal_type<T, std::string>::value
				|| boost::is_arithmetic<T>::value>::type load(T& value,
		boost::property_tree::ptree const& pt, string const& key,
		bool complete) {
	if (complete || pt.find(key) != pt.not_found())
		load(value, pt.get_child(key), complete);
}

template<class T>
void load_items(std::vector<T>& vec, boost::property_tree::ptree const& pt,
		string const& key, bool complete) {
	string vector_key = key + string(".count");
	if (complete || pt.find(vector_key) != pt.not_found()) {
		size_t count = pt.get<size_t>(vector_key);

		for (size_t i = 0; i != count; ++i) {
			T t;
			load(t, pt.get_child(str(format("%s.item_%d") % key % i)),
					complete);
			vec.push_back(t);
		}
	}
}

void inline split(vector<string>& vec, string const& space_separated_list) {
	std::istringstream iss(space_separated_list);
	while (iss) {
		std::string value;
		iss >> value;
		if (value.length()) {
			vec.push_back(value);
		}
	}
}

void inline load_split(vector<string>& vec, boost::property_tree::ptree const& pt, string const& key) {
	boost::optional<string> values = pt.get_optional<string>(key);
	if (values) {
		split(vec, *values);
	}
}

template<class T>
void inline load(vector<T>& vec, boost::property_tree::ptree const& pt, string const& key, bool complete) {
	boost::optional<T> value = pt.get_optional<T>(key);
	if (value) {
		vec.push_back(*value);
		return;
	}
	for (size_t i = 1;; i++) {
		value = pt.get_optional<string>(key + "#" + ToString(i));
		if (value) {
			vec.push_back(*value);
			continue;
		}
		value = pt.get_optional<string>(key + "." + ToString(i));
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
void load(T& value, boost::property_tree::ptree const& pt, string const& key) {
	load(value, pt, key, true);
}

template<class T>
void load(T& value, boost::property_tree::ptree const& pt, const char* key) {
	load(value, pt, string(key), true);
}

template<class T>
void load(T& value, boost::property_tree::ptree const& pt) {
	load(value, pt, true);
}

// config singleton-wrap
template<class Config>
struct config {
    static std::string dirnameOf(const std::string& fname)
    {
        size_t pos = fname.find_last_of("\\/");
        return (std::string::npos == pos) ? "" : fname.substr(0, pos);
    }


	static void create_instance(std::string const& filename) {
		boost::property_tree::ptree pt;
		boost::property_tree::read_info(filename, pt);
        load(inner_cfg(), pt);
		is_initialized() = true;
	}

	static Config const& get() {
		VERIFY(is_initialized());
		return inner_cfg();
	}

	static Config& get_writable() {
		VERIFY(is_initialized());
		return inner_cfg();
	}

private:
	static Config& inner_cfg() {
		static Config config;
		return config;
	}

	static bool& is_initialized() {
		static bool is_initialized = false;
		return is_initialized;
	}
};

}

template<class T>
inline void load_param(const string& filename, const string& key,
		boost::optional<T>& value) {
	boost::property_tree::ptree pt;
	boost::property_tree::read_info(filename, pt);
	boost::optional<T> loaded_value = pt.get_optional<T>(key);
	value = loaded_value;
}

template<class T>
inline void write_param(const string& filename, const string& key,
		const optional<T>& value) {
	if (value) {
		std::ofstream params_stream(filename.c_str(), std::ios_base::app);
		params_stream << key << "\t" << value << std::endl;
	}
}

template<class K, class V>
inline void load_param_map(const string& filename, const string& key,
		map<K, V>& value) {
	boost::property_tree::ptree pt;
	boost::property_tree::read_info(filename, pt);
	boost::optional<std::string> as_str = pt.get_optional<std::string>(key);
	if (as_str) {
		vector<std::string> key_value_pairs;
		boost::split(key_value_pairs, *as_str, boost::is_any_of(";"));
		for (auto it = key_value_pairs.begin(); it != key_value_pairs.end();
				++it) {
			vector<std::string> key_value;
			boost::split(key_value, *it, boost::is_any_of(" "));
			VERIFY(key_value.size() == 2);
			value[boost::lexical_cast<K>(key_value[0])] =
					boost::lexical_cast<K>(key_value[1]);
		}
	}
}

template<class K, class V>
inline void write_param_map(const string& filename, const string& key,
		const map<K, V>& value) {
	if (value.size() > 0) {
		ofstream params_stream(filename.c_str(), std::ios_base::app);
		params_stream << key << "\t\"";
		string delim = "";
		for (auto it = value.begin(); it != value.end(); ++it) {
			params_stream << delim << it->first << " " << it->second;
			delim = ";";
		}
		params_stream << "\"" << endl;
	}
}

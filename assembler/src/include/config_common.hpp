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
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "simple_tools.hpp"
#include "fs_path_utils.hpp"

namespace config_common
{
    // for enable_if/disable_if
	namespace details
	{
		template <class T, class S>
		struct is_equal_type
		{
			static const bool value = false;
		};

		template <class T>
		struct is_equal_type<T, T>
		{
			static const bool value = true;
		};
	}

	void correct_relative_includes(const fs::path& p, fs::path const& working_dir = fs::initial_path())
    {
	    using namespace boost;
	    using namespace boost::filesystem;

	    std::ifstream f(p.string().c_str());

	    vector<string> strings;

	    while (f.good())
	    {
	        string raw_line;
	        getline(f, raw_line);

	        strings.push_back(raw_line);
	        boost::trim(raw_line);

	        const string incl_prefix = "#include";

	        if (starts_with(raw_line, incl_prefix))
	        {
	            string incl_path_str = raw_line.substr(incl_prefix.size(), string::npos);
	            trim   (incl_path_str);
	            trim_if(incl_path_str, is_any_of("\""));

	            path incl_path = p.parent_path() / incl_path_str;
	                 incl_path = make_relative_path(incl_path);

	            correct_relative_includes(incl_path);

	            strings.back() = incl_prefix + " \"" + incl_path.string() + "\"";
	        }
	    }

	    f.close();

	    std::ofstream of (p.string().c_str());
	    for (size_t i = 0; i < strings.size(); ++i)
	        of << strings[i] << std::endl;
    }


    template <class T>
    typename boost::enable_if_c<details::is_equal_type<T, std::string>::value || boost::is_arithmetic<T>::value >::type
        load(T& value, boost::property_tree::ptree const& pt, string const& key, bool complete)
    {
        if (complete || pt.find(key) != pt.not_found())
            value = pt.get<T>(key);
    }

    template <class T>
    typename boost::disable_if_c<details::is_equal_type<T, std::string>::value || boost::is_arithmetic<T>::value >::type
        load(T& value, boost::property_tree::ptree const& pt, string const& key, bool complete)
    {
        if (complete || pt.find(key) != pt.not_found())
            load(value, pt.get_child(key), complete);
    }

	template<class T>
	void load(std::vector<T>& vec, boost::property_tree::ptree const& pt, string const& key, bool complete)
	{
		string vector_key = key + string(".count");
		if (complete || pt.find(vector_key) != pt.not_found())
		{
			size_t count = pt.get < size_t > (vector_key);

			for (size_t i = 0; i != count; ++i)
			{
				T t;
				load(t, pt.get_child(str(format("%s.item_%d") % key % i)), complete);
				vec.push_back(t);
			}
		}
	}

    template<class T>
    void load(T& value, boost::property_tree::ptree const& pt, string const& key)
    {
        load(value, pt, key, true);
    }

    template<class T>
    void load(T& value, boost::property_tree::ptree const& pt, const char* key)
    {
        load(value, pt, string(key), true);
    }

    template<class T>
    void load(T& value, boost::property_tree::ptree const& pt)
    {
        load(value, pt, true);
    }


//    template<class T>
//    void load(T&, boost::property_tree::ptree const&, bool complete);

	// config singleton-wrap
	template <class Config>
	struct config
	{
	//	template<typenamea...Args>
	//	static void create_instance(std::string const& filename, Args&&... args)
	//	{
	//		boost::property_tree::ptree pt;
	//		boost::property_tree::read_info(filename, pt);
	//		load(pt, inner_cfg(), std::forward<Args>(args)...);
	//	}

		static void create_instance(std::string const& filename)
		{
			correct_relative_includes(filename);

			boost::property_tree::ptree pt;
			boost::property_tree::read_info(filename, pt);
			load(inner_cfg(), pt);
		}

		static Config const& get()
		{
			return inner_cfg();
		}

		static Config& get_writeable()
		{
			return inner_cfg();
		}

	private:
		static Config& inner_cfg()
		{
			static Config config;
			return config;
		}
	};

}

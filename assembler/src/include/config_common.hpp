/*
 * config_common.hpp
 *
 *  Created on: Aug 13, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_COMMON_HPP_
#define CONFIG_COMMON_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <string>
#include <vector>
#include <iostream>
#include "simple_tools.hpp"

using namespace std;
using boost::format;
//using boost::str;

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
				load(t, pt.get_child(str(format("%$s.item_%d") % key % i)));
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
			boost::property_tree::ptree pt;
			boost::property_tree::read_info(filename, pt);
			load(inner_cfg(), pt);
		}

		static Config const& get()
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

#endif /* CONFIG_COMMON_HPP_ */

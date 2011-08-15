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
#include <string>

// for enable_if/disable_if
namespace
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
typename boost::enable_if_c<is_equal_type<T, std::string>::value || boost::is_arithmetic<T>::value >::type
	load(boost::property_tree::ptree const& pt, std::string const& key, T& value)
{
	value = pt.get<T>(key);
}

template <class T>
typename boost::disable_if_c<is_equal_type<T, std::string>::value || boost::is_arithmetic<T>::value >::type
	load(boost::property_tree::ptree const& pt, std::string const& key, T& value)
{
	load(pt.get_child(key), value);
}

template<class T>
void load(boost::property_tree::ptree const&, T&);


// config singleton-wrap
template <class Config>
struct config
{
	static void create_instance(std::string const& filename)
	{
		boost::property_tree::ptree pt;
		boost::property_tree::read_info(filename, pt);
		load(pt, inner_cfg());
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


#endif /* CONFIG_COMMON_HPP_ */

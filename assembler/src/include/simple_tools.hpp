//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * simple_tools.hpp
 *
 *  Created on: 27.05.2011
 *      Author: vyahhi
 */

#ifndef SIMPLE_TOOLS_HPP_
#define SIMPLE_TOOLS_HPP_

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "verify.hpp"
#include "io/ireader.hpp"
#include <memory>
#include <string>
#include <set>
#include <vector>

#include <boost/format.hpp>

/**
 * Converts anything to string (using ostringstream).
 */
template <typename T>
std::string ToString(const T& t) {
	std::ostringstream ss;
	ss << t;
	return ss.str();
}

template <typename T>
std::string ToString(const T& t, size_t length) {
	std::ostringstream ss;
	ss << t;
	std::string result = ss.str();
	while(result.size() < length)
		result = "0" + result;
	return result;
}

template <typename T>
std::string ToString(std::vector<T>& t) {
	std::ostringstream ss;
	ss << "Size "<<t.size()<<": [";
	for (auto it = t.begin(); it != t.end(); ++it)
		ss<<*it<<", ";
	ss<<"]";
	return ss.str();
}

template <typename T>
std::string ToString(std::set<T>& t) {
	std::ostringstream ss;
	ss << "Size "<<t.size()<<": [";
	for (auto it = t.begin(); it != t.end(); ++it)
		ss<<*it<<", ";
	ss<<"]";
	return ss.str();
}

//taken from http://habrahabr.ru/post/131977/
class FormattedString {

 public:
  FormattedString(const char* fmt): m_fmt(fmt)
  {
  }

  template<class T>
  FormattedString& operator<< (const T& arg) {
    m_fmt % arg;
    return *this;
  }

  operator std::string() const {
    return m_fmt.str();
  }

  std::string str() const {
    return m_fmt.str();
  }

 protected:
  boost::format m_fmt;
};

template<typename T>
inline const std::pair<T, T> ReversePair(std::pair<T, T> ep) {
  return std::pair<T, T>(ep.second, ep.first);
}

template <class ContainerT1, class ContainerT2>
void push_back_all(ContainerT1& target, const ContainerT2& to_insert) {
	target.insert(target.end(), to_insert.begin(), to_insert.end());
}

template <class ContainerT1, class ContainerT2>
void insert_all(ContainerT1& target, const ContainerT2& to_insert) {
	target.insert(to_insert.begin(), to_insert.end());
}

template<class MapT>
std::set<typename MapT::key_type> key_set(const MapT& m) {
    std::set<typename MapT::key_type> answer;
    for (auto it = m.begin(); it != m.end(); ++it) {
        answer.insert(it->first);
    }
    return answer;
}

template<class MapT>
std::set<typename MapT::mapped_type> value_set(const MapT& m) {
    std::set<typename MapT::mapped_type> answer;
    for (auto it = m.begin(); it != m.end(); ++it) {
        answer.insert(it->second);
    }
    return answer;
}

template <class MapT>
const typename MapT::mapped_type& get(const MapT& from, const typename MapT::key_type& key) {
	auto it = from.find(key);
	VERIFY(it != from.end());
	return it->second;
}

template <class MapT>
typename MapT::mapped_type& get(MapT& from, const typename MapT::key_type& key) {
	auto it = from.find(key);
	VERIFY(it != from.end());
	return it->second;
}

template <class MMapT>
const std::vector<typename MMapT::mapped_type> get_all(const MMapT& from, const typename MMapT::key_type& key) {
    std::vector<typename MMapT::mapped_type> answer;
    for (auto it = from.lower_bound(key); it != from.upper_bound(key); ++it) {
        answer.push_back(it->second);
    }
	return answer;
}

namespace std
{
template<class T1, class T2>
std::ostream& operator<< (std::ostream& os, std::pair<T1, T2> const& pair)
{
	return os << "(" << pair.first << ", " << pair.second << ")";
}
//}

//namespace omnigraph
//{
template<class T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v)
{
 	os << "[";
 	std::string delim = "";
 	for (auto it = v.begin(); it != v.end(); ++it) {
 		os << delim << *it;
 		delim = ", ";
 	}
// 	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, ", "));
 	os << "]";
 	return os;
}

}

#endif /* SIMPLE_TOOLS_HPP_ */

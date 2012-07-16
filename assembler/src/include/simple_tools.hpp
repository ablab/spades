//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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

#include "verify.hpp"
#include "io/ireader.hpp"
#include <memory>
#include <string>
#include <set>
#include <vector>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

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

//template <typename T>
//std::string str(const T& t) {
//	return ToString(t);
//}
//
//template <typename T>
//std::string str(std::vector<T>& t) {
//	return ToString(t);
//}

template<class T>
std::auto_ptr<T> make_ptr(T* t) {
	return std::auto_ptr<T>(t);
}

template<class T>
std::auto_ptr<T> make_shared_ptr(T* t) {
	return boost::shared_ptr<T>(t);
}

//taken from http://habrahabr.ru/post/131977/
class FormattedString {
public:
	FormattedString(const char* fmt): m_fmt(fmt) {}

    template<class T>
    FormattedString& operator<< (const T& arg) {
        m_fmt % arg;
        return *this;
    }
    operator std::string() const {
        return m_fmt.str();
    }
protected:
    boost::format m_fmt;
};

/**
 * Checks if file exists.
 * Analogs: http://www.techbytes.ca/techbyte103.html , http://www.gamedev.net/topic/211918-determining-if-a-file-exists-c/
 */
bool fileExists(std::string filename);

inline bool fileExists(std::string filename) {
	namespace fs = boost::filesystem;
	return fs::is_regular_file(filename);
}

/**
 * Exit(1) if file doesn't exists, writes FATAL log message.
 */
inline void checkFileExistenceFATAL(std::string filename) {
	if (!fileExists(filename)) {
		VERIFY_MSG(false, "File " << filename << " doesn't exist or can't be read!\n");
	}
}

template <class ContainerT1, class ContainerT2>
void insert_all(ContainerT1& target, const ContainerT2& to_insert) {
	target.insert(target.end(), to_insert.begin(), to_insert.end());
}

template <class map_t>
const typename map_t::_Tp& get(const map_t& from, const typename map_t::_Key& key) {
	auto it = from.find(key);
	VERIFY(it != from.end());
	return it->second;
}

namespace std
{
template<class T1, class T2>
std::ostream& operator<< (std::ostream& os, std::pair<T1, T2> const& pair)
{
	return os << "(" << pair.first << ", " << pair.second << ")";
}
}

namespace omnigraph
{
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

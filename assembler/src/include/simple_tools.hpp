/*
 * simple_tools.hpp
 *
 *  Created on: 27.05.2011
 *      Author: vyahhi
 */

#ifndef SIMPLE_TOOLS_HPP_
#define SIMPLE_TOOLS_HPP_

#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include "logging.hpp"
#include "verify.hpp"
#include "io/ireader.hpp"

#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem.hpp>

#include <fstream>

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

template<class T>
std::auto_ptr<T> create_auto_ptr(T* t) {
	return std::auto_ptr<T>(t);
}

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

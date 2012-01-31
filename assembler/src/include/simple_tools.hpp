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

/**
 * Checks if file exists.
 * Analogs: http://www.techbytes.ca/techbyte103.html , http://www.gamedev.net/topic/211918-determining-if-a-file-exists-c/
 */
bool fileExists(std::string filename);

/**
 * Exit(1) if file doesn't exists, writes FATAL log message.
 */
void checkFileExistenceFATAL(std::string filename);

/**
 * Use vector<T> as input-stream with operator>>(T& t)
 */
template <typename T>
class VectorStream : public io::IReader<T> {
	std::vector<T> data_;
	size_t pos_;
	bool closed_;
public:
	VectorStream(const std::vector<T>& data) : data_(data), pos_(0), closed_(false) {

	}

	virtual ~VectorStream() {

	}

	virtual bool eof() /*const */{
		return pos_ == data_.size();
	}

	VectorStream<T>& operator>>(T& t) {
		t = data_[pos_++];
		return *this;
	}

	void close() {
		closed_ = true;
	}

	virtual bool is_open() /*const */{
		return !closed_;
	}

	void reset() {
		pos_ = 0;
	}

};

inline bool fileExists(std::string filename) {
	namespace fs = boost::filesystem;
	return fs::is_regular_file(filename);
}

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
std::ostream& operator<< (std::ostream& os, std::vector<T> const& v)
{
 	os << "[";
 	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, ", "));
 	os << "]";
 	return os;
}
}

#endif /* SIMPLE_TOOLS_HPP_ */

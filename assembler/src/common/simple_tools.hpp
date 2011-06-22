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
#include <vector>

/**
 * Converts anything to string (using ostringstream).
 */
template <typename T>
std::string ToString(T& t) {
	std::ostringstream ss;
	ss << t;
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
class VectorStream {
	std::vector<T> data_;
	size_t pos_;
	bool closed_;
public:
	VectorStream(const std::vector<T>& data) : data_(data), pos_(0), closed_(false) {

	}

	bool eof() const {
		return pos_ == data_.size();
	}

	VectorStream<T>& operator>>(T& t) {
		t = data_[pos_++];
		return *this;
	}

	void close() {
		closed_ = true;
	}

	bool is_open() const {
		return !closed_;
	}

	void reset() {
		pos_ = 0;
	}

};

#endif /* SIMPLE_TOOLS_HPP_ */

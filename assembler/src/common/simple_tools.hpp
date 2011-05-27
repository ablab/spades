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
 * Use vector<T> as input-stream with operator>>(T& t)
 */
template <typename T>
class VectorStream {
	vector<T> data_;
	size_t pos_;
	bool closed_;
public:
	VectorStream(const vector<T>& data) : data_(data), pos_(0), closed_(false) {

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

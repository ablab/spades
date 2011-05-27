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

template <typename T>
std::string ToString(T t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

template <typename T>
class VectorStream {
	vector<T> data_;
	size_t pos_;
	bool closed_;
public:
	VectorStream(const vector<T>& data) : data_(data), pos_(0), closed_(false) {

	}

	bool eof() {
		return pos_ == data_.size();
	}

	VectorStream<T>& operator>>(T& t) {
		t = data_[pos_++];
		return *this;
	}

	void close() {
		closed_ = true;
	}

	bool is_open() {
		return !closed_;
	}

	void reset() {
		pos_ = 0;
	}

};

#endif /* SIMPLE_TOOLS_HPP_ */

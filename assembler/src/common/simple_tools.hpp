#ifndef SIMPLE_TOOLS_HPP_
#define SIMPLE_TOOLS_HPP_

#include <sstream>
#include <string>
#include <vector>

//namespace common {

using std::string;
using std::vector;

string Reverse(const string &s);

string Complement(const string &s);

string ReverseComplement(const string &s);

template <typename T>
string ToString(T t) {
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
//}
#endif /* SIMPLE_TOOLS_HPP_ */

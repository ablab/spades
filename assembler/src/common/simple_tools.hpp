#ifndef SIMPLE_TOOLS_HPP_
#define SIMPLE_TOOLS_HPP_

#include <string>
using std::string;

string Reverse(const string &s);

string Complement(const string &s);

string ReverseComplement(const string &s);

template <typename T>
string toString(T t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

#endif /* SIMPLE_TOOLS_HPP_ */

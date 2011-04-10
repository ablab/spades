/*
 * test_utils.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

template<typename T>
struct PairHash {
	size_t operator()(pair<T, T> p) const {
		return hash<T> ()(p.first) + hash<T> ()(p.second);
	}
};

template<typename T>
struct PairLess {
	bool operator()(pair<T, T> p1, pair<T, T> p2) const {
		return less<T> ()(p1.first, p2.first) ? true : (less<T> ()(p2.first,
				p1.first) ? false : less<T> ()(p1.second, p2.second));
	}
};

std::string complement(const std::string& s) {
	return (!Sequence(s)).str();
}

vector<Read> MakeReads(string *ss, size_t count) {
	vector<Read> ans;
	for (size_t i = 0; i < count; ++i) {
		Read r("", *ss, "");
		ss++;
		ans.push_back(r);
	}
	return ans;
}

#endif /* TEST_UTILS_HPP_ */

/*
 * strobe_read.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef STROBE_READ_HPP_
#define STROBE_READ_HPP_

#include "seq.hpp"
#include <iostream>
using namespace std;

/*
 * Use single_read<size,T>::type for single read
 * Use mate_read<size,T>::type for mate read
 * Use strobe_read<size,cnt,T>::type for cnt-read
 *
 * where size -- number of nucleotides in each read
 * (we don't store info about gaps' size here)
 */
template <int size, int cnt = 1, typename T = char>
class strobe_read {
public:
	strobe_read() {
		// random sequences constructor
	}

	void put(int i, const string &s) { // by value (probably faster)
		data_[i] = Seq<size,T>(s);
	}

	Seq<size,T> operator[](size_t i) {
		return data_[i];
	}

private:
	array<Seq<size,T>,cnt> data_;
};

// use mate_read<size, T>::type for mate reads
template <int size, typename T = char>
struct mate_read {
	typedef strobe_read<size,2,T> type; // because "There is no direct way to have templated typedefs in C++" :(
};

// use single_read<size, T>::type for single reads
template <int size, typename T = char>
struct single_read {
	typedef strobe_read<size,1,T> type; // because "There is no direct way to have templated typedefs in C++" :(
};

#endif /* STROBE_READ_HPP_ */

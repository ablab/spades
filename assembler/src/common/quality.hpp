/*
 * qual.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef QUAL_HPP_
#define QUAL_HPP_

#include <array>
#include <string>
using namespace std;

class Quality {
public:

	Quality(const string &s) : qual_(s) {
		//
	}

	int operator[](size_t i) const {
		return qual_[i];
	}

	string str() const { // copying (defensive)!
		return qual_;
	}

private:
	string qual_;
	friend class ireadstream;
};

#endif /* QUAL_HPP_ */

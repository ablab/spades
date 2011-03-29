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
	Quality(const string &s) : qual(s) {
	}
	int operator[](size_t i) {
		return qual[i];
	}
private:
	string qual;
};

#endif /* QUAL_HPP_ */

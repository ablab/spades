/*
 * qual.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef QUAL_HPP_
#define QUAL_HPP_

template <int size>
class Qual {
public:
	Qual(const string &s) {
		for (size_t i = 0; i < size; ++i) {
			quals[i] = s[i];
		}
	}
	int operator[](size_t i) {
		return quals[i];
	}
private:
	std::array<char,size> quals;
};

#endif /* QUAL_HPP_ */

/*
 * Simple operations and checks for nucleotide-letters
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#ifndef NUCL_HPP_
#define NUCL_HPP_

#include <cassert>

inline bool is_nucl(char c) { // is ACGT
	return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
}

inline bool is_dignucl(char c) { // is 0123
	return (c >= 0 && c < 4);
}

inline char complement(char c) { // 0123 -> 3210
	assert(is_dignucl(c));
	return c ^ 3;
}

inline char nucl(char c) { // 0123 -> ACGT
	assert(is_dignucl(c));
	switch(c) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N'; // never happens
	}
}

inline char dignucl(char c) { // ACGT -> 0123
	assert(is_nucl(c));
	switch(c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return -1; // never happens
	}
}


#endif /* NUCL_HPP_ */

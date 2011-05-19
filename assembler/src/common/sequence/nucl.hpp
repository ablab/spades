/**
 * @file    nucl.hpp
 * @author  vyahhi
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Simple operations and checks for nucleotide-letters
 *
 */


#ifndef NUCL_HPP_
#define NUCL_HPP_

#include <cassert>

/**
* @param char c
* @return true if c is 'A', 'C', 'G', 'T'.
*/
inline bool is_nucl(char c) { // is ACGT
	return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
}

inline bool is_dignucl(char c) { // is 0123
	return (c >= 0 && c < 4);
}

/**
* 0123 -> 3210
* @param char c
* @return c ^ 3
*/
inline char complement(char c) {
	assert(is_dignucl(c));
	return c ^ 3;
}

/**
* @param char c
* @return complement symbol
*/
inline char nucl_complement(char c){
	assert(is_nucl(c));
	switch(c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default: return -1; // never happens
	}
}

/**
* 0123 -> ACGT
* @param char c
* @return complement symbol
*/
inline char nucl(char c) {
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
		default: assert(false);//return -1; // never happens
	}
}


#endif /* NUCL_HPP_ */

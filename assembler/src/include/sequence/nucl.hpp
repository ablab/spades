//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

#include "verify.hpp"
#include <iostream>

// some optimization: instead of switches we use predefined constant arrays... surprisingly works faster (10% on input).

const char dignucl_map['T' + 1] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3};

const bool isnucl_map[256] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

const char nucl_map[4] = {'A', 'C', 'G', 'T'};

const char nucl_complement_map['T' + 1] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A'};

/**
 * ACGT -> true
 * @param char c
 * @return true if c is 'A', 'C', 'G' or 'T'.
 */
inline bool is_nucl(signed char c) { // is ACGT
	return isnucl_map[c];
}

/**
 * 0123 -> true
 * @param char c
 * @return true if c is 0, 1, 2 or 3.
 */
inline bool is_dignucl(char c) { // is 0123
	return (c < 4 && c >= 0);
}

/**
 * 0123 -> 3210
 * @param char c
 * @return c ^ 3
 */
inline char complement(char c) {
	//VERIFY(is_dignucl(c));
	return c ^ 3;
}

/**
 * ACGT -> TGCA
 * @param char c is 'A', 'C', 'G' or 'T'
 * @return complement symbol, i.e. 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'.
 */

struct nucl_complement_functor { // still unused
	inline bool operator() (signed char c) const {
		return nucl_complement_map[c];
	}
};

inline char nucl_complement(signed char c){
	//VERIFY(is_nucl(c));
	return nucl_complement_map[c];
	/*
	switch(c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		//case 'N': return 'N';
		default: VERIFY(false); return -1; // never happens
	}*/
}

/**
 * 0123 -> ACGT
 * @param char c is 0, 1, 2 or 3
 * @return 0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T'
 */
inline char nucl(signed char c) {
	return nucl_map[c];
	/*
	VERIFY(is_dignucl(c));
	switch(c) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N'; // never happens
	}
	*/
}

/**
 * ACGT -> 0123
 * @param char c is 'A', 'C', 'G' or 'T'
 * @return A => 0, C => 1, G => 2, T => 3
 */

/*
struct dignucl : public unary_function<int,bool> {
	bool operator()(signed char c) const {
		return dignucl_map[c];		
	}	
};*/

inline char dignucl(signed char c) {
	return dignucl_map[c];
	//VERIFY(is_nucl(c));
	/*switch(c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: VERIFY(false); return -1; // never happens
	}*/
}


#endif /* NUCL_HPP_ */

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

#include "utils/verify.hpp"
#include <iostream>

const char dignucl_map['T' + 1] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3};

const bool isnucl_map[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

const char nucl_map[4] = {'A', 'C', 'G', 'T'};

const  char nucl_complement_map['T' + 1] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A'};

/**
 * ACGT -> true
 * @param char c
 * @return true if c is 'A', 'C', 'G' or 'T'.
 */
inline bool is_nucl(char c) { // is ACGT
    return isnucl_map[(unsigned)c];
}

/**
 * 0123 -> true
 * @param char c
 * @return true if c is 0, 1, 2 or 3.
 */
inline bool is_dignucl(char c) { // is 0123
    return (c < 4);
}

/**
 * 0123 -> 3210
 * @param char c
 * @return c ^ 3
 */
inline char complement(char c) {
    // VERIFY(is_dignucl(c));
    return c ^ 3;
}

/**
 * ACGT -> TGCA
 * @param char c is 'A', 'C', 'G', 'T' or 'N'
 * @return complement symbol, i.e. 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N'
 */

struct nucl_complement_functor { // still unused
    inline bool operator() (char c) const {
        char cc = nucl_complement_map[(unsigned)c];
        return cc ? cc : 'N';
    }
};

inline char nucl_complement(char c){
    // TODO: deal with 'N' case
    //VERIFY(is_nucl(c));
    char cc = nucl_complement_map[(unsigned)c];
    return cc ? cc : 'N';
}

/**
 * 0123 -> ACGT
 * @param char c is 0, 1, 2 or 3
 * @return 0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T'
 */
inline char nucl(char c) {
    return nucl_map[(unsigned)c];
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

inline char dignucl(char c) {
    // VERIFY(is_nucl(c));
    return dignucl_map[(unsigned)c];
}


#endif /* NUCL_HPP_ */

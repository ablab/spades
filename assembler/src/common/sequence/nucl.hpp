//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef NUCL_HPP_
#define NUCL_HPP_

#include "utils/verify.hpp"

/**
 * 0123 -> true
 * @param char c
 * @return true if c is 0, 1, 2 or 3.
 */
inline bool is_dignucl(char c) {
    return (c < 4);
}

/**
 * 0123 -> 3210
 * @param char c
 * @return c ^ 3
 */
inline char complement(char c) {
    VERIFY_DEV(is_dignucl(c));
    return char(c ^ 3);
}

static const char INVALID_NUCL = char(-1);

/**
 * ACGTacgt0123 -> true
 * @param char c
 * @return true if c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3'.
 */

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

// WARNING: These functions were carefully crafted for speed
// Change only if you know what you're doing
inline bool is_nucl(char c) {
    if (unlikely(c>= 0 && c < 4))
        return true;

    switch (c) {
        case 'a':
        case 'A':
        case 'c':
        case 'C':
        case 'g':
        case 'G':
        case 't':
        case 'T':
            return true;
        default:
            return false;
    }
}

/**
 * ACGT -> TGCA
 * @param char c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3' or 'N'
 * @return complement symbol, i.e. 'A/a/0' => 'T/t/3', 'C/c/1' => 'G/g/2', 'G/g/2' => 'C/c/1', 'T/t/3' => 'A/a/0', 'N' => 'N'
 */
inline char nucl_complement(char c) {
    if (unlikely(c>= 0 && c < 4))
        return complement(c);

    switch (c) {
        case 'a':
            return 't';
        case 'A':
            return 'T';
        case 'c':
            return 'g';
        case 'C':
            return 'G';
        case 'g':
            return 'c';
        case 'G':
            return 'C';
        case 't':
            return 'a';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        case 'n':
            return 'n';
        default:
            VERIFY_DEV(false);
            return INVALID_NUCL;
    }
}

/**
 * 0123acgtACGT -> ACGT
 * @param char c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3'
 * @return 'A/a/0' => 'A', 'C/c/1' => 'C', 'G/g/2' => 'G', 'T/t/3' => 'T'
 */
inline char nucl(char c) {
    if (likely(c >= 0 && c < 4)) {
        const uint32_t nucl_map = 'A' + ('C' << 8) + ('G' << 16) + ('T' << 24);
        return (char)((nucl_map >> (8*c)) & 0xFF);
    } else if ('A' <= c && c <= 'T')
        return c;

    return (char)(c - 'a' + 'A');
}

/**
 * ACGT -> 0123
 * @param char c is 'A/a/', 'C', 'G' or 'T'
 * @return A => 0, C => 1, G => 2, T => 3
 */
inline char dignucl(char c) {
    if (unlikely(c>= 0 && c < 4)) {
        return c;
    } else if (unlikely('a' <= c && c <= 't')) {
        c = (char)(c - 'a' + 'A');
    }

    return (c <= 'C' ?
            (c == 'A' ? 0 : 1) :
            (c == 'G' ? 2 : 3));
}

#undef likely
#undef unlikely

#endif /* NUCL_HPP_ */

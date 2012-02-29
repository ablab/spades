// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#ifndef SHUFFLE_H
#define SHUFFLE_H

#include "Vec.h"



/**
 * Shuffle
 *
 * Given an integer N>0, it will fill shuffled with the shuffled integers
 * between 0 and N-1 (included).
 *
 * It uses std::random_shuffle and a function object which uses
 * drand48_r. Thus, it is thread safe and multiprocessing safe, in the
 * sense that the same seed will always produce the same sequence.
 */
void Shuffle( int N, vec<int> &shuffled, int seed = 0 );
void Shuffle64( uint64_t N, vec<uint64_t> &shuffled, uint64_t seed = 0 );


#endif

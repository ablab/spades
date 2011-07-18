/**
 * @file    kmer_functions.hpp
 * @author  adavydow, snikolenko
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
 * In this there are several functions which can be helpfull while
 * generating kmers from reads.
 *
 */
#ifndef HAMMER_KMERFUNCTIONS_HPP_
#define HAMMER_KMERFUNCTIONS_HPP_
#include <vector>
#include "common/read/read.hpp"
#include "hammer/kmer_stat.hpp"
#warning("This file is deprecated. For fast k-mer iteration use valid_kmer_generator.hpp")
/**
 * trim bad quality nucleotides from start and end of the read
 * @return size of the read left
 */
uint32_t TrimBadQuality(Read *r, int bad_quality_threshold = 2) __attribute__ ((deprecated));

/**
 * Deprecated: Will be replaced soon.
 * @param k k as in k-mer
 * @param start start point
 * @return the first starting point of a valid k-mer >=start; return
 * seq_.size() if no such place exists
 */
uint32_t FirstValidKmerPos(const Read &r, uint32_t start, uint32_t k) __attribute__ ((deprecated));

Sequence GetSubSequence(const Read &r, uint32_t start, uint32_t length) __attribute__ ((deprecated));

template<uint32_t kK>
vector< Seq<kK> > GetKMers(const Read &r) __attribute__ ((deprecated));

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v) __attribute__ ((deprecated));

/**
 * Deprecated: Syntax is to hard and bug-dangerous. While such a
 * construction can be rather useful for performance it's better to
 * write a class providing same functionality.
 * @param kmer get next valid k-mer
 * @param pos starting point
 * @return the first starting point of a valid k-mer >=start; return
 * -1 if no such place exists
 */
template<uint32_t kK>
int NextValidKmer(const Read &r, int prev_pos, Seq<kK> *kmer) __attribute__ ((deprecated));

#include "kmer_functions.impl.hpp" //  Template function implementation
#endif  // HAMMER_KMERFUNCTIONS_HPP_

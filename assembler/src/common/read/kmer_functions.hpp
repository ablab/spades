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
#include "hammer/hammer_config.hpp"

class Read;
class Sequence;
/**
 * trim bad quality nucleotides from start and end of the read
 * @return size of the read left
 */
uint32_t TrimBadQuality(Read *r, int bad_quality_threshold = 2);

/**
 * @param k k as in k-mer
 * @param start start point
 * @return the first starting point of a valid k-mer >=start; return
 * seq_.size() if no such place exists
 */
uint32_t FirstValidKmerPos(const Read &r, uint32_t start, uint32_t k);

/**
 * @param kmer get next valid k-mer
 * @param pos starting point
 * @return the first starting point of a valid k-mer >=start; return
 * -1 if no such place exists
 */
int32_t NextValidKmer(const Read &r, int32_t pos, KMer & kmer);

/**
 * add k-mers from read to map
 */
void AddKMers(const Read &r, KMerStatMap *v);

/**
  * Finds next valid kmer in a read.
  * @return position of the next valid kmer in a read; returns -1 if no more valid kmers
  * @param pos should equal -1 if it's the first time and start pos of kmer otherwise
  */
size_t NextValidKmer(const Read &r, size_t pos, KMer & kmer);

Sequence GetSubSequence(const Read &r, uint32_t start, uint32_t length);
#endif  // HAMMER_KMERFUNCTIONS_HPP_

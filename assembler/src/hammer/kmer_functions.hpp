/*
* kmer_functions.hpp
*
*  Created on: 13.07.2011
*      Author: adavydow
*/
#ifndef HAMMER_KMERFUNCTIONS_HPP_
#define HAMMER_KMERFUNCTIONS_HPP_
#include "hammer/hammer_config.hpp"
#include "common/read/read.hpp"

/**
 * trim bad quality nucleotides from start and end of the read
 * @return size of the read left
 */
size_t TrimBadQuality(Read &r, int bad_quality_threshold = 2);

/**
 * @param k k as in k-mer
 * @param start start point
 * @return the first starting point of a valid k-mer >=start; return seq_.size() if no such place exists
 */
size_t FirstValidKmerPos(const Read &r, size_t start, size_t k);

/**
 * add k-mers from read to map
 */
void AddKMers(const Read &r, KMerStatMap &v);

Sequence GetSubSequence(const Read &r, size_t start, size_t length);
#endif // HAMMER_KMERFUNCTIONS_HPP_

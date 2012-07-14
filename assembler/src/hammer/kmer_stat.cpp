//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_stat.cpp
 *
 *  Created on: 10.03.2012
 *      Author: snikolenko
 */

#include "standard.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"

char getQual(const KMerCount & kmc, size_t i) {
	if (Globals::use_common_quality)
    return Globals::common_quality * kmc.second.count;
	if (kmc.second.count == 1)
		return Globals::blobquality[kmc.first.start() + i];
	else
    return kmc.second.qual[i];
}

double getProb(const KMerCount &kmc, size_t i) {
  uint8_t qual = getQual(kmc, i);

  return Globals::quality_probs[qual];
}

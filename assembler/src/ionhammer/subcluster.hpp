//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __SUBCLUSTER_HPP__
#define __SUBCLUSTER_HPP__

#include <vector>
#include <cstddef>

#include "hkmer.hpp"

class KMerData;

size_t subcluster(KMerData &kmer_data, std::vector<size_t> &cluster);

// for debug purposes
hammer::HKMer center(const KMerData &data, const std::vector<size_t>& kmers);

#endif // __SUBCLUSTER_HPP__

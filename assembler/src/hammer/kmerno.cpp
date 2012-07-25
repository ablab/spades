//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "position_kmer.hpp"
#include "kmerno.hpp"

bool KMerNo::operator==(const KMerCount &kmc) const {
  return ((strncmp(Globals::blob + info.index, Globals::blob + kmc.first.start(), K) == 0));
}

std::string KMerNo::str() const {
  return std::string(Globals::blob + info.index, K);
}

uint64_t KMerNo::hash::operator()(const KMerNo &kn) const {
  size_t h = 239;
  for (size_t i = 0; i < K; i++) {
    h = ((h << 5) - h) + Globals::blob[kn.info.index+i];
  }
  return h;
}

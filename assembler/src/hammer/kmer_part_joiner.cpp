/**
 * @file    kmer_part_joiner.cpp
 * @author  Alex Davydow
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
 * This class provides functionality which allows us to join
 * several sorted files containing k-mer and its frequency to the one
 * file, still sorted.
 */
#include "kmer_part_joiner.hpp"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

using std::pair;
using std::string;
using std::swap;
using std::vector;

KMerPartJoiner::KMerPartJoiner(const vector<FILE*> &ifiles, int k)
    : kmer_parsers_(), k_(k) {
  for (size_t i = 0; i < ifiles.size(); ++i) {
    KMerPartParser kpp(ifiles[i], k_);
    if (!kpp.eof()) {
      kmer_parsers_.insert(kpp);
    }
  }
}

pair<string, int> KMerPartJoiner::Next() {
  KMerPartParser kpp(*kmer_parsers_.begin());
  pair<string, int> ret = make_pair(kpp.last_string(), kpp.last_count());
  kmer_parsers_.erase(kpp);
  kpp.Next();
  if (!kpp.eof()) {
    kmer_parsers_.insert(kpp);
  }
  return ret;
}

void KMerPartJoiner::KMerPartParser::Swap(KMerPartParser other) {
  swap(file_, other.file_);
  swap(last_string_, other.last_string_);
  swap(last_count_, other.last_count_);
  swap(eof_, other.eof_);
}

void KMerPartJoiner::KMerPartParser::Next() {
  char *buf = new char[k_ + 1];
  char format[10];
  snprintf(format, sizeof(format), "%%%ds %%d", k_);
  eof_ = (fscanf(file_, format, buf, &last_count_) == EOF);
  if (!eof_) {
    last_string_ = buf;
  }
  delete[] buf;
}

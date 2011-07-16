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
#include <utility>
#include <string>
#include <vector>

using std::vector;
using std::pair;
using std::string;

KMerPartJoiner::KMerPartJoiner(const vector<FILE*> &ifiles) {
  for (size_t i = 0; i < ifiles.size(); ++i) {
    KMerPartParser kpp(ifiles[i]);
    if (!kpp.eof()) {
      kmer_parsers_.insert(kpp);
      }
  }
}

pair<string, int> KMerPartJoiner ::Next() {
  KMerPartParser kpp(*kmer_parsers_.begin());
  pair<string, int> ret = make_pair(kpp.last_string(), kpp.last_count());
  kmer_parsers_.erase(kpp);
  kpp.Next();
  if (!kpp.eof()) {
    kmer_parsers_.insert(kpp);
  }
  return ret;
}

bool KMerPartJoiner::IsEmpty() {
  return kmer_parsers_.size() == 0;
}

KMerPartJoiner::KMerPartParser::KMerPartParser(FILE *file) {
  file_ = file;
  eof_ = false;
  Next();
}

bool KMerPartJoiner::KMerPartParser::operator<(
                                     const KMerPartParser &other) const {
  return last_string_ < other.last_string_;
}

KMerPartJoiner::KMerPartParser::KMerPartParser(const KMerPartParser &other) {
  file_ = other.file_;
  last_string_ = other.last_string_;
  last_count_ = other.last_count_;
  eof_ = other.eof_;
}

void KMerPartJoiner::KMerPartParser::Next() {
  char buf[K + 1];
  eof_ = (fscanf(file_, "%s %d", buf, &last_count_) == EOF);
  last_string_ = buf;
}

bool KMerPartJoiner::KMerPartParser::eof() {
  return eof_;
}

string KMerPartJoiner::KMerPartParser::last_string() {
  return last_string_;
}

int KMerPartJoiner::KMerPartParser::last_count() {
  return last_count_;
}

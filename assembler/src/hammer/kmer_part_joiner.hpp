/**
 * @file    kmer_part_joiner.hpp
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
#ifndef HAMMER_KMERPARTJOINER_HPP_
#define HAMMER_KMERPARTJOINER_HPP_

#include <utility>
#include <vector>
#include <set>
#include <string>
#include <cstdio>
#include "hammer/hammer_config.hpp"

class KMerPartJoiner {
 public:
  explicit KMerPartJoiner(const std::vector<FILE*> &ifiles);
  std::pair<string, int> Next();
  bool IsEmpty();
 private:

  class KMerPartParser {
   public:
    explicit KMerPartParser(FILE *file);
    bool operator<(const KMerPartParser &other) const;
    KMerPartParser(const KMerPartParser &other);
    void Next();
    bool eof();
    std::string last_string();
    int last_count();
   private:
    std::string last_string_;
    int last_count_;
    FILE *file_;
    bool eof_;
  };

  std::set<KMerPartParser> kmer_parsers_;
};

#endif  // HAMMER_KMERPARTJOINER_HPP_

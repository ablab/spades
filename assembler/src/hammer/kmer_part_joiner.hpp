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
 * This class provides functionality which allows us to join several 
 * sorted files containing k-mer and its frequency to the one file, 
 * still sorted.
 */
#ifndef HAMMER_KMERPARTJOINER_HPP_
#define HAMMER_KMERPARTJOINER_HPP_

#include <cstdio>
#include <utility>
#include <vector>
#include <set>
#include <string>
#include "hammer/kmer_freq_info.hpp"

class KMerPartJoiner {
 public:
  explicit KMerPartJoiner(const std::vector<FILE*> &ifiles, int k);
  KMerFreqInfo Next();
  bool IsEmpty() const { return kmer_parsers_.size() == 0; }
 private:
  class KMerPartParser {
   public:
    KMerPartParser(const KMerPartParser &other)
        : last_(),
          file_(other.file_),
          eof_(other.eof_),
          k_(other.k_) {}
    KMerPartParser& operator=(const KMerPartJoiner::KMerPartParser other) {
      Swap(other);
      return *this;
    }
    explicit KMerPartParser(FILE *file, int k)
        : last_(),
          file_(file),
          eof_(false),
          k_(k){
      Next();
    }
    bool eof() const { return eof_; }
    KMerFreqInfo last() {
      return last_;
    }
    void Next();
    class Lesser {
    public:
      bool operator()(const KMerPartParser &a, const KMerPartParser &b) const {
        return a.last_.kmer < b.last_.kmer;
      }
    };
   private:
    void Swap(KMerPartParser other);
    KMerFreqInfo last_;
    FILE *file_;
    bool eof_;
    int k_;
  };
  std::set<KMerPartParser, KMerPartParser::Lesser> kmer_parsers_;
  int k_;
  // Disallow copy and assign.
  KMerPartJoiner(const KMerPartJoiner&);
  void operator=(const KMerPartJoiner&);
};

#endif  // HAMMER_KMERPARTJOINER_HPP_

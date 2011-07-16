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
#include "hammer/hammer_config.hpp"

class KMerPartJoiner {
 public:
  explicit KMerPartJoiner(const std::vector<FILE*> &ifiles);
  std::pair<string, int> Next();
  bool IsEmpty() const { return kmer_parsers_.size() == 0; }
 private:
  class KMerPartParser {
   public:
    KMerPartParser(const KMerPartParser &other)
        : last_string_(other.last_string_),
          last_count_(other.last_count_),
          file_(other.file_),
          eof_(other.eof_) {}
    KMerPartParser& operator=(const KMerPartJoiner::KMerPartParser other) {
      Swap(other);
      return *this;
    }
    explicit KMerPartParser(FILE *file)
        : last_string_(""),
          last_count_(-1),
          file_(file),
          eof_(false) {
      Next();
    }
    bool eof() const { return eof_; }
    std::string last_string() const { return last_string_; }
    int last_count() const { return last_count_; }
    void Next();
    class Lesser {
    public:
      bool operator()(const KMerPartParser &a, const KMerPartParser &b) const {
        return a.last_string_ < b.last_string_;
      }
    };
   private:
    void Swap(KMerPartParser other);
    std::string last_string_;
    int last_count_;
    FILE *file_;
    bool eof_;
  };
  std::set<KMerPartParser, KMerPartParser::Lesser> kmer_parsers_;
  // Disallow copy and assign.
  KMerPartJoiner(const KMerPartJoiner&);
  void operator=(const KMerPartJoiner&);
};

#endif  // HAMMER_KMERPARTJOINER_HPP_

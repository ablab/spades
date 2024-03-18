
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#undef NDEBUG
#include <algorithm>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
#include <locale>

// to compile: gcc this_prog.c -lz
#include <stdio.h>
#include <zlib.h>

#include <kseq/kseq.h>
KSEQ_INIT(gzFile, gzread)

#include "fasta_reader.hpp"
#include "utils/logger/logger.hpp"
#include "utils/stl_utils.hpp"

std::string rev_comp(const std::string &s) {
  std::string result;
  std::unordered_map<char, char> rc = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}, {'N', 'N'}};
  result.reserve(s.length());
  for (auto it = s.rbegin(), end = s.rend(); it != end; ++it) {
    auto map_it = rc.find(*it);
    if (map_it != rc.cend()) {
      result += map_it->second;
    } else {
      result += 'N';
      WARN("Unrecognized symbol: " << *it);
    }
  }
  return result;
}

std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc) {
  std::vector<std::string> edges;

  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(filename.c_str(), "r");
  assert(fp);
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    // printf("name: %s\n", seq->name.s);
    // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    // printf("seq: %s\n", seq->seq.s);
    // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
    std::string read = seq->seq.s;

    // U -> T
    for (char &ch : read) {
      if (ch == 'U') {
        ch = 'T';
      }
    }

    edges.push_back(read);
    if (add_rc) {
      edges.push_back(rev_comp(read));
    }
  }

  kseq_destroy(seq);
  gzclose(fp);

  return edges;
}

std::vector<std::pair<std::string, std::string>> read_fasta(const std::string &filename) {
  std::vector<std::pair<std::string, std::string>> records;

  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(filename.c_str(), "r");
  assert(fp);
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string name(seq->name.s);
    std::string s(seq->seq.s);
    utils::trim(name);
    utils::trim(s);
    records.emplace_back(std::move(name), std::move(s));
  }

  kseq_destroy(seq);
  gzclose(fp);

  return records;
}

// vim: set ts=2 sw=2 et :

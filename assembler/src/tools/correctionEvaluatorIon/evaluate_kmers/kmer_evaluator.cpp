//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "hkmer.hpp"
#include "io/reader.hpp"
#include "io/single_read.hpp"
#include "io/read_processor.hpp"

#include <cstddef>
#include <iostream>
#include <mutex>
#include <unordered_set>

void printUsage() {
  std::cerr << "usage: ./kmer_evaluator <reference.fasta> <contigs.fasta>" << std::endl;
}

namespace std {
  template<>
  struct hash<hammer::HSeq<16ul> > {
    size_t operator()(hammer::HSeq<16ul> seq) const {
      return seq.GetHash();
    }
  };
}

typedef std::unordered_set<hammer::HKMer> HKMerSet;

struct SetFiller {
  HKMerSet& kmers;
  std::mutex mutex_;

  SetFiller(HKMerSet& kmers) : kmers(kmers) {}

  bool operator()(const io::SingleRead& read) {
    processString(read.GetSequenceString());
    return false;
  }

  void processString(const std::string &seq) {
    if (seq.empty())
      return;
    std::vector<hammer::HKMer> kmers;
    kmers.reserve(seq.size());
    hammer::HKMer kmer;
    char nucl = seq[0];
    kmer <<= nucl;
    size_t n_runs = 1;
    for (size_t i = 1; i < seq.size(); ++i) {
      if (seq[i] != nucl) {
        if (n_runs++ >= hammer::K) {
          kmers.push_back(kmer);
          kmers.push_back(!kmer);
        }
        nucl = seq[i];
      }
      kmer <<= seq[i];
    }
    kmers.push_back(kmer);
    kmers.push_back(!kmer);
    pushKMers(kmers);
  }

  void pushKMers(const std::vector<hammer::HKMer> &hkmers) {
    mutex_.lock();
    for (auto it = hkmers.begin(); it != hkmers.end(); ++it) {
      auto& hkmer = *it;
      kmers.insert(hkmer);
    }
    mutex_.unlock();
  }
};

void fillSet(HKMerSet& kmers, char* filename) {
  const unsigned num_threads = 16;

  SetFiller filler(kmers);
  io::Reader irs(filename, io::PhredOffset);
  hammer::ReadProcessor(num_threads).Run(irs, filler);
}

void runComparison(const HKMerSet& reference_kmers, const HKMerSet& contig_kmers) {
  size_t total_genomic = reference_kmers.size();
  size_t total_contig = contig_kmers.size();

  size_t contig_genomic = 0;

  for (auto it = reference_kmers.cbegin(), et = reference_kmers.cend(); it != et; ++it)
    if (contig_kmers.find(*it) != contig_kmers.end())
      ++contig_genomic;
  
  long contig_non_genomic = total_contig - contig_genomic;

  std::cout << "Reference kmers:    " << total_genomic << std::endl;
  std::cout << "Contig kmers:       " << total_contig << std::endl;
  std::cout << "  Genomic:          " << contig_genomic 
            << " (" << (contig_genomic * 100.0 / total_genomic) << "%)" << std::endl;
  std::cout << "  Non-genomic:      " << contig_non_genomic << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    printUsage(); return 0;
  }

  HKMerSet reference, contigs;
  reference.reserve(10000000);
  contigs.reserve(200000000);
  std::cout << "Filling set of reference kmers..." << std::endl;
  fillSet(reference, argv[1]);
  std::cout << "Filling set of contig kmers..." << std::endl;
  fillSet(contigs, argv[2]);
  runComparison(reference, contigs);
}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <memory>
#include "hkmer.hpp"
#include "io/reads/read_processor.hpp"
#include "kmer_helpers.h"

void printUsage() {
  std::cerr << "usage: ./kmer_evaluator <reference.fasta> <contigs.fasta>"
            << std::endl;
}

void runComparison(const HKMerSet& reference_kmers,
                   const HKMerSet& contig_kmers) {
  size_t total_genomic = reference_kmers.size();
  size_t total_contig = contig_kmers.size();

  size_t contig_genomic = 0;

  for (auto it = reference_kmers.cbegin(), et = reference_kmers.cend();
       it != et; ++it)
    if (contig_kmers.find(*it) != contig_kmers.end()) ++contig_genomic;

  long contig_non_genomic = total_contig - contig_genomic;

  std::cout << "Reference kmers:    " << total_genomic << std::endl;
  std::cout << "Contig kmers:       " << total_contig << std::endl;
  std::cout << "  Genomic:          " << contig_genomic << " ("
            << ((double)contig_genomic * 100.0 / (double)total_genomic) << "%)" << std::endl;
  std::cout << "  Non-genomic:      " << contig_non_genomic << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    printUsage();
    return 0;
  }

  HKMerSet reference, contigs;
  reference.reserve(10000000);
  contigs.reserve(200000000);
  std::cout << "Filling set of reference kmers..." << std::endl;
  FillSet(reference, argv[1]);
  std::cout << "Filling set of contig kmers..." << std::endl;
  FillSet(contigs, argv[2]);
  std::cout << "Running comparison " << std::endl;
  runComparison(reference, contigs);
}

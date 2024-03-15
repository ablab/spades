//
// Created by Vasiliy Ershov on 10/07/16.
//

#ifndef PROJECT_REFERENCE_H
#define PROJECT_REFERENCE_H
#include "utils/logger/log_writers.hpp"

#include "hkmer.hpp"

#include "io/reads/file_reader.hpp"
#include "io/reads/read_processor.hpp"
#include "kmer_helpers.h"

#include <cstddef>
#include <iostream>
#include <mutex>
#include <unordered_set>

class TGenomReferenceOracle {
 private:
  const std::string FilePath;
  HKMerSet ReferenceKMers;

 public:
  TGenomReferenceOracle(const std::string& filePath) : FilePath(filePath) {
    FillSet(ReferenceKMers, filePath.data());
    INFO("Reference kmers:    " << ReferenceKMers.size());
  }

  bool IsGenomic(const hammer::HKMer& kmer) const {
    return ReferenceKMers.count(kmer) > 0;
  }

  void KMerSetStats(const HKMerSet& kmers, std::string setName) const {
    INFO("Stats for " << setName);

    size_t total_genomic = ReferenceKMers.size();
    size_t total_set = kmers.size();

    size_t set_genomic = 0;

    for (auto it = ReferenceKMers.cbegin(), et = ReferenceKMers.cend();
         it != et; ++it) {
      if (kmers.count(*it) > 0) {
        set_genomic += 1;
      }
    }

    long set_non_genomic = total_set - set_genomic;

    INFO("Set kmers:       " << total_set);
    INFO("Genomic: " << set_genomic << " ("
                     << ((double)set_genomic * 100.0 / (double)total_genomic) << "%)");
    INFO("NonGenomic: " << set_non_genomic);
  }
};

#endif  // PROJECT_REFERENCE_H

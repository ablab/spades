//
// Created by Vasiliy Ershov on 10/07/16.
//

#ifndef PROJECT_KMER_HELPERS_H
#define PROJECT_KMER_HELPERS_H

#include <mutex>
#include <unordered_set>
#include "hkmer.hpp"
#include "io/reads/file_reader.hpp"
#include "io/reads/read_processor.hpp"
#include "valid_hkmer_generator.hpp"

using HKMerSet = std::unordered_set<hammer::HKMer>;

namespace std {
template <>
struct hash<hammer::HSeq<hammer::K> > {
  size_t operator()(hammer::HSeq<hammer::K> seq) const { return seq.GetHash(); }
};
}  // namespace std

class SetFiller {
 private:
  std::unordered_set<hammer::HKMer>& kmers_;
  std::mutex mutex_;

 private:
  void ProcessString(const std::string& seq) {
    if (seq.empty()) {
      return;
    }
    std::vector<hammer::HKMer> kmers;
    kmers.reserve(seq.size());
    ValidHKMerGenerator<hammer::K> generator(seq.data(), nullptr, seq.size());
    while (generator.HasMore()) {
      kmers.push_back(generator.kmer());
      kmers.push_back(!generator.kmer());
      generator.Next();
    }
    PushKMers(kmers);
  }

  void PushKMers(const std::vector<hammer::HKMer>& hkmers) {
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto it = hkmers.begin(); it != hkmers.end(); ++it) {
      auto& hkmer = *it;
      kmers_.insert(hkmer);
    }
  }

 public:
  SetFiller(std::unordered_set<hammer::HKMer>& kmers) : kmers_(kmers) {}

  bool operator()(std::unique_ptr<io::SingleRead>&& read) {
    ProcessString(read->GetSequenceString());
    return false;
  }
};

inline void FillSet(HKMerSet& kmers, const char* filename) {
  const unsigned num_threads = 16;
  SetFiller filler(kmers);
  io::FileReadStream irs(filename, io::PhredOffset);
  hammer::ReadProcessor(num_threads).Run(irs, filler);
}

#endif  // PROJECT_KMER_HELPERS_H

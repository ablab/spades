//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef CCLEAN_ADAPTERINDEX_HPP
#define CCLEAN_ADAPTERINDEX_HPP

#include "sequence/seq.hpp"
#include "utils/mph_index/kmer_index.hpp"

#include <string>
#include <set>
#include <unordered_map>

namespace cclean {
const unsigned K = 10;
typedef Seq<K> KMer;

class AdapterIndex {
  typedef std::set<std::size_t> IndexValueType;
  std::unordered_map<KMer, IndexValueType, KMer::hash> index_;

 public:
  AdapterIndex() {}

  void clear() {
    index_.clear();
    seqs_.clear();
  }
  IndexValueType& operator[](cclean::KMer s) { return index_[s]; }
  auto find(cclean::KMer s) const -> decltype(index_.find(s)) { return index_.find(s); }
  auto end() const -> decltype(index_.end()) { return index_.end(); }

  bool contains(cclean::KMer s) const {
    return index_.find(s) != index_.end();
  }
  const std::string& seq(size_t idx) const { return seqs_[idx]; }

 private:
  std::vector<std::string> seqs_;

  friend class AdapterIndexBuilder;
};

class AdapterIndexBuilder {
 public:
  AdapterIndexBuilder() {}

  void FillAdapterIndex(const std::string &db, AdapterIndex &index);

 private:
  DECL_LOGGER("Index Building");
};

  // end of namespace
}

#endif // __CCLEAN__ADAPTERINDEX_HPP__

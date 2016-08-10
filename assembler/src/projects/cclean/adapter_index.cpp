//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "adapter_index.hpp"
#include "io/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/ireadstream.hpp"
#include "config_struct_cclean.hpp"

#include <libcxx/sort.hpp>

using namespace cclean;

void AdapterIndexBuilder::FillAdapterIndex(const std::string &db, AdapterIndex &data) {
  data.clear();

  INFO("Reading adapter database from " << db);
  ireadstream irs(db);
  while (!irs.eof()) {
    Read r;
    irs >> r;
    const std::string &seq = r.getSequenceString();

    data.seqs_.push_back(seq);
    data.seqs_.push_back(ReverseComplement(seq));
  }

  INFO("Filling adapter index");
  for (size_t i = 0, e = data.seqs_.size(); i !=e; ++i) {
    const std::string &seq = data.seqs_[i];
    ValidKMerGenerator<cclean::K> gen(seq.c_str(), NULL, seq.size());

    while (gen.HasMore()) {
      KMer kmer = gen.kmer();

      auto& entry = data.index_[kmer];
      entry.insert(i);

      gen.Next();
    }
  }

  INFO("Done. Total " << data.seqs_.size() << " adapters processed. Total "
                      << data.index_.size() << " unique k-mers.");
}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_IO_READ_CORRECTOR_HPP__
#define __HAMMER_IT_IO_READ_CORRECTOR_HPP__

#include "HSeq.hpp"
#include "config_struct.hpp"
#include "consensus.hpp"
#include "flow_space_read.hpp"
#include "hkmer_distance.hpp"
#include "valid_hkmer_generator.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/optional.hpp>

#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/SamHeader.h>
#include "kmer_data.hpp"
#include "seqeval/BaseHypothesisEvaluator.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <vector>

#if 1
#include <iomanip>
#include <iostream>
#endif

namespace hammer {
namespace correction {

template <class ReadCorrector>
class SingleReadCorrector {
  const KMerData& hkmer_data_;
  using PenaltyCalcer = typename ReadCorrector::PenaltyCalcer;
  using Factory = typename PenaltyCalcer::PenaltyCalcerFactory;

 public:
  struct ReadSelectionPredicate {
    virtual bool operator()(const io::SingleRead& read) = 0;
  };

  struct DebugOutputPredicate : public ReadSelectionPredicate {};

  struct NoDebug : public DebugOutputPredicate {
    virtual bool operator()(const io::SingleRead&) { return false; }
  };

  struct FullDebug : public DebugOutputPredicate {
    virtual bool operator()(const io::SingleRead&) { return true; }
  };

  class DebugIfContains : public DebugOutputPredicate {
    Sequence Needle;
    Sequence NeedleRc;

   public:
    DebugIfContains(const Sequence& seq) : Needle(seq), NeedleRc(!seq) {}

    virtual bool operator()(const io::SingleRead& read) {
      auto readSeq = read.sequence();
      if (readSeq.size() < Needle.size()) return false;
      if (readSeq.find(Needle, 0) != -1ULL) return true;
      return readSeq.find(NeedleRc, 0) != -1ULL ? true : false;
    }
  };

  struct SelectPredicate : public ReadSelectionPredicate {};

  struct SelectAll : public SelectPredicate {
    virtual bool operator()(const io::SingleRead&) { return true; }
  };

  class SelectByName : public SelectPredicate {
    std::set<std::string> Names;

   public:
    SelectByName(const std::set<std::string>& names) : Names(names) {}

    virtual bool operator()(const io::SingleRead& r) {
      return Names.find(r.name()) != Names.end();
    }
  };

 private:
  BamTools::SamHeader* sam_header_ptr_;
  DebugOutputPredicate& debug_predicate_;
  SelectPredicate& select_predicate_;
  const Factory& penalty_factory_;
  ReadCorrector read_corrector_;

 public:
  SingleReadCorrector(const KMerData& kmer_data,
                      const Factory& penalty_factory,
                      BamTools::SamHeader* sam_header,
                      DebugOutputPredicate& debug, SelectPredicate& select)
      : hkmer_data_(kmer_data),
        sam_header_ptr_(sam_header),
        debug_predicate_(debug),
        select_predicate_(select),
        penalty_factory_(penalty_factory),
        read_corrector_(hkmer_data_, penalty_factory_) {}

  SingleReadCorrector(const KMerData& kmer_data,
                      const Factory& penalty_factory,
                      DebugOutputPredicate& debug, SelectPredicate& select)
      : hkmer_data_(kmer_data),
        sam_header_ptr_(NULL),
        debug_predicate_(debug),
        select_predicate_(select),
        penalty_factory_(penalty_factory),
        read_corrector_(hkmer_data_, penalty_factory_) {}

  std::unique_ptr<io::SingleRead> operator()(std::unique_ptr<io::SingleRead> r) {
    return SingleReadCorrector::operator()(*r);
  }

  std::unique_ptr<io::SingleRead> operator()(const io::SingleRead& read) {
    if (!select_predicate_(read)) {
      return nullptr;
    }

    bool debug_mode = debug_predicate_(read);
    if (debug_mode) {
      std::cerr << "=============================================" << std::endl;
      std::cerr << '>' << read.name() << '\n'
                << read.GetSequenceString() << std::endl;
    }
    auto corected_seq = read_corrector_.Correct(
        read, cfg::get().keep_uncorrected_ends, debug_mode);

    if (corected_seq.empty()) {
      return nullptr;
    }

    auto result = std::unique_ptr<io::SingleRead>(new io::SingleRead(read.name(), corected_seq));
    return result;
  }

  std::unique_ptr<io::BamRead> operator()(std::unique_ptr<BamTools::BamAlignment> alignment) {
    VERIFY(sam_header_ptr_);
    io::SingleRead r(alignment->Name, alignment->QueryBases);
    // reverse strand means we're working with a mapped BAM, might be
    // the case for datasets downloaded from IonCommunity
    if (alignment->IsReverseStrand()) r = !r;
    auto corrected_r = SingleReadCorrector::operator()(r);
    std::string rg;
    if (!alignment->GetTag("RG", rg) || !corrected_r) return nullptr;
    auto flow_order = sam_header_ptr_->ReadGroups[rg].FlowOrder;

    float delta_score, fit_score;
    auto seq = corrected_r->GetSequenceString();
    if (alignment->IsReverseStrand()) {
      std::reverse(seq.begin(), seq.end());
      for (auto it = seq.begin(); it != seq.end(); ++it) {
        switch (*it) {
          case 'A':
            *it = 'T';
            break;
          case 'C':
            *it = 'G';
            break;
          case 'G':
            *it = 'C';
            break;
          case 'T':
            *it = 'A';
            break;
          default:
            break;
        }
      }
    }

    BaseHypothesisEvaluator(*alignment, flow_order, seq, delta_score, fit_score,
                            0);
    std::stringstream ss;
    ss << alignment->Name << "_" << delta_score << "_" << fit_score;
    alignment->Name = ss.str();
    if (delta_score >= cfg::get().delta_score_threshold)
      return std::unique_ptr<io::BamRead>(new io::BamRead(*alignment));

    BamTools::BamAlignment corrected(*alignment);
    corrected.QueryBases = corrected_r->GetSequenceString();
    return std::unique_ptr<io::BamRead>(new io::BamRead(corrected));
  }
};

template <class ReadCorrector>
class PairedReadCorrector : public SingleReadCorrector<ReadCorrector> {
public:

  using PenaltyCalcer = typename ReadCorrector::PenaltyCalcer;
  using Factory = typename PenaltyCalcer::PenaltyCalcerFactory;

 public:
  PairedReadCorrector(
      const KMerData& kmerData, const Factory& penaltyFactory,
      typename SingleReadCorrector<ReadCorrector>::DebugOutputPredicate& debug,
      typename SingleReadCorrector<ReadCorrector>::SelectPredicate& select)
      : SingleReadCorrector<ReadCorrector>(kmerData, penaltyFactory, debug,
                                            select) {}

  std::unique_ptr<io::PairedRead> operator()(std::unique_ptr<io::PairedRead> r) {
    auto corrected_r = SingleReadCorrector<ReadCorrector>::operator()(r->first());
    auto corrected_l = SingleReadCorrector<ReadCorrector>::operator()(r->second());

    if (!corrected_r || !corrected_l) return nullptr;

    return std::unique_ptr<io::PairedRead>(
        new io::PairedRead(*corrected_r, *corrected_l, 0));
  }
};

};      // namespace correction
};      // namespace hammer
#endif  // __HAMMER_IT_IO_READ_CORRECTOR_HPP__

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

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

#include "read_corrector_structs_new.h"

namespace hammer {
namespace correction {

template <class CorrectionsLikelihoodCalcer>
class ReadCorrector {
 public:
  using PenaltyCalcer = CorrectionsLikelihoodCalcer;
 private:
  using State = CorrectionState<typename PenaltyCalcer::PenaltyState>;
  const KMerData& data;
  using PenaltyCalcerFactory = typename CorrectionsLikelihoodCalcer::PenaltyCalcerFactory;
  const PenaltyCalcerFactory& penalty_calcer_factory;

  mutable size_t skipped_reads = 0;
  mutable size_t queue_overflow_reads = 0;

  inline bool Flush(std::priority_queue<State>& candidates,
                    std::priority_queue<State>& corrections,
                    size_t limit,
                    size_t readSize) const {

    if (corrections.size() > limit) {
      auto top = pop_queue(candidates);
      if (!std::isinf(top.Penalty())) {
        corrections.emplace(std::move(top));
      }
      std::priority_queue<State>().swap(candidates);
      return true;
    } else {
      while (!candidates.empty()) {
        auto top = pop_queue(candidates);
        if (top.TotalCorrections() > std::max(readSize / 10, (size_t)3)) {
          continue;
        }
        if (!std::isinf(top.Penalty())) {
          corrections.emplace(std::move(top));
        }
      }
      return false;
    }
  }

  std::string CorrectRight(const PenaltyCalcer& penalty_calcer,
                           const std::string& read,
                           const size_t offset,
                           bool reverse,
                           bool& is_too_many_corrections,
                           bool make_only_simple_corrections = false) const {
    if (offset >= read.size()) {
      return read;
    }

    std::priority_queue<State> corrections;
    std::priority_queue<State> candidates;

    CorrectionContext context(data, read, reverse);
    {
      corrections.emplace(StateBuilder<PenaltyCalcer>::Initial(
          context, penalty_calcer, (uint)offset));
    }

    std::map<uint, std::set<size_t> > visited;
    const size_t queue_limit =  (const size_t)(cfg::get().queue_limit_multiplier * log2(read.size() - offset + 1));//(const size_t)(100 * read.size());

    bool queue_overflow = false;

    while (!corrections.empty()) {

      auto state = std::pop_queue(corrections);
      assert(state.Position() <= read.size());

      {
        size_t hash = state.GetHKMer().GetHash();
        if (visited[state.Position()].count(hash) && corrections.size()) {
          continue;
        }
        visited[state.Position()].insert(hash);
      }

      if (state.Position() < read.size()) {
        MoveToNextDivergence<PenaltyCalcer> mover(state,
                                                  context,
                                                  penalty_calcer);
        if (mover.FindNextDivergence()) {
          mover.Move();
        }
      }

      if (state.Position() == read.size()) {
        return state.Read()->ToString();
      }

      //      //don't correct last kmer
      if ((state.Position() + context.GetHRun(state.Position()).len) ==
          read.size()) {
        auto result = state.Read()->ToString();
        result += (context.GetHRun(state.Position()).str());
        return result;
      }

      {
        SkipMayBeBadHRun<PenaltyCalcer> skipHRun(state,
                                                 context,
                                                 penalty_calcer);
        candidates.emplace(skipHRun.State());
      }

      {
        CorrectLastHRun<PenaltyCalcer> hrun_corrector(state,
                                                      context,
                                                      penalty_calcer);
        if (make_only_simple_corrections) {
          hrun_corrector.AddOnlySimpleCorrections(candidates);
        } else {
          hrun_corrector.AddPossibleCorrections(candidates);
        }
        queue_overflow |= Flush(candidates, corrections, queue_limit, read.size());
      }
    }
    is_too_many_corrections = queue_overflow;

    return read;
  }

 public:

  ReadCorrector(const KMerData& kmer_data,
                 const PenaltyCalcerFactory& factory)
      : data(kmer_data)
      , penalty_calcer_factory(factory) {}

  ~ReadCorrector() {
    INFO("Skipped reads count: " << skipped_reads);
    if (queue_overflow_reads) {
      WARN("Too many possible corrections in some reads (" << queue_overflow_reads << "), something may be wrong");
    }
  }

  std::string Correct(const io::SingleRead& read,
                      bool keep_uncorrected_ends = true,
                      bool debug = false,
                      uint simple_passes_count = 0,
                      uint complex_passes_count = 1) const {

    std::string current_read = read.GetSequenceString();

    PenaltyCalcer penalty_calcer = penalty_calcer_factory(current_read);

    bool overflow = false;

    for (uint pass = 0; pass < 2 * (simple_passes_count + complex_passes_count); ++pass) {
      const bool reverse = pass % 2 == 0;  // tail has more errors, so let's start with "simple" part
      const bool only_simple = pass < 2 * simple_passes_count;
      if (reverse) {
        current_read = ReverseComplement(current_read);
      }
      const auto solid_island = penalty_calcer.SolidIsland(current_read);
      const size_t solid_length = solid_island.right_ - solid_island.left_;

      if (debug) {
#pragma omp critical
        {
          std::cerr << "Solid length: " << solid_length << " / "
                    << current_read.size() << std::endl;
          std::cerr << "Position: " << solid_island.left_ << " / "
                    << solid_island.right_ << std::endl;
        }
      }

      if (solid_length == 0 || solid_length == current_read.size()) {
        if (pass == 0) {
          if (solid_length == 0) {
#pragma omp atomic
            skipped_reads++;
          }
        }

        break;
      }

      bool pass_overflow = false;
      current_read = CorrectRight(penalty_calcer,
                                  current_read,
                                  solid_island.right_,
                                  reverse,
                                  overflow,
                                  only_simple);

      overflow  |= pass_overflow;

      if (reverse) {
        current_read = ReverseComplement(current_read);
      }
    }

    if (overflow) {
        #pragma omp atomic
        queue_overflow_reads++;
    }

    if (!keep_uncorrected_ends) {
      return penalty_calcer.TrimBadQuality(current_read);
    }
    return current_read;
  }
};

};      // namespace correction
};      // namespace hammer
#endif  // __HAMMER_IT_READ_CORRECTOR_HPP__

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef ADDITIONAL_CPP
#define ADDITIONAL_CPP

#include "output.hpp"
#include "config_struct_cclean.hpp"
#include "io/read_processor.hpp"

  enum WorkModeType {
    NONE = 0,
    SINGLE_END = 1,
    SINGLE_END_Q = 2,
    BRUTE_SIMPLE = 3,
    BRUTE_WITH_Q = 4
  };

  constexpr double MatchScore = 0.6;
  constexpr double MismatchScore = 100;

  class AbstractCclean {
      // Abstract base class for cclean functors
    public:
      AbstractCclean(std::ostream &aligned_output, std::ostream &bed,
                     const std::string &db,
                     const WorkModeType &mode,
                     const unsigned mlen,
                     const bool full_inform = false)
                :aligned_(0), full_inform_(full_inform), read_mlen_(mlen),
                 mismatch_threshold_(cfg::get().mismatch_threshold),
                 score_threshold_(cfg::get().score_treshold),
                 aligned_part_fraction_(cfg::get().aligned_part_fraction),
                 db_name_(db), mode_(mode), aligned_output_stream_(aligned_output),
                 bad_stream_(bed)  {}
      virtual Read operator()(const Read &read, bool *ok) = 0;
      inline size_t aligned() { return aligned_; }
      virtual ~AbstractCclean() {}

    protected:
      size_t aligned_;

      const bool full_inform_;
      const uint read_mlen_;
      const double mismatch_threshold_;  // for nonquality mode
      const double score_threshold_;  // for quality mode

      const double aligned_part_fraction_;
      const std::string &db_name_;
      const WorkModeType mode_;

      std::ostream &aligned_output_stream_;
      std::ostream &bad_stream_;
      // Abstract for clean functors
      class AbstractCleanFunctor {
        public:
          inline virtual bool operator()(const Read &r,
                          const StripedSmithWaterman::Alignment &a,
                          double aligned_part, const std::string &adapter,
                          double *best_score) = 0;
          virtual ~AbstractCleanFunctor() {}
      };
  };

#endif // ADDITIONAL_CPP

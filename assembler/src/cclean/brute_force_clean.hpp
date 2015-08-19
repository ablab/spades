//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef BRUTE_FORCE_CLEAN_HPP
#define BRUTE_FORCE_CLEAN_HPP

#include "utils.hpp"
#include "additional.cpp"

class BruteForceClean: public AbstractCclean {
  // Class that get read with oper() and clean it, if that possible
  public:
    BruteForceClean(std::ostream& aligned_output,
                    std::ostream& bed,const std::string &db,
                    const WorkModeType &mode,
                    const uint mlen,
                    const std::vector<std::string> &gen,
                    const bool full_inform = false)
      : AbstractCclean(aligned_output, bed, db, mode, mlen, full_inform),
        adap_seqs_(gen)  {
      if(mode == BRUTE_SIMPLE) checker = new BruteCleanFunctor;
      if(mode == BRUTE_WITH_Q) checker = new BruteQualityCleanFunctor;
    }
    virtual ~BruteForceClean() { delete checker; }
    // ReadProcessor class put each read in this operator
    virtual Read operator()(const Read &read, bool *ok);

  private:
    const std::vector<std::string> &adap_seqs_;
    std::string best_adapter_;
    AbstractCleanFunctor *checker; // Checks is adapter in read

    // Here goes functors for clean in different modes
    class BruteCleanFunctor: public AbstractCleanFunctor {
        virtual inline bool operator()(const Read &r,
                                       const StripedSmithWaterman::Alignment &a,
                                       double aligned_part, const std::string &adapter,
                                       double *best_score) {
          double cur_score = cclean_utils::
                             GetMismatches(r.getSequenceString(), adapter, a);
          if (cur_score < (*best_score) &&
              cclean_utils::is_alignment_good(a, r.getSequenceString(), adapter,
                                              aligned_part)) {
            (*best_score) = cur_score;
            return true;
          }
          return false;
        }
    };
    class BruteQualityCleanFunctor: public AbstractCleanFunctor {
        virtual inline bool operator()(const Read &r,
                                       const StripedSmithWaterman::Alignment &a,
                                       double aligned_part, const std::string &adapter,
                                       double *best_score) {
          double cur_score = cclean_utils::
                             GetScoreWithQuality(a, r.getQuality().str());
          if (cur_score >= (*best_score) &&
              cclean_utils::is_alignment_good(a, r.getSequenceString(), adapter,
                                              aligned_part)) {
            (*best_score) = cur_score;
            return true;
          }
          return false;
        }
    };
};

#endif // BRUTE_FORCE_CLEAN_HPP

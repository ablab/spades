//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef JOB_WRAPERS_HPP
#define JOB_WRAPERS_HPP

#include "additional.cpp"
#include "utils.hpp"

namespace cclean {
  class AdapterIndex;
}

class SimpleClean: public AbstractCclean {
  public:
    SimpleClean(std::ostream &aligned_output,
                std::ostream &bed, const std::string &db,
                const WorkModeType &mode,
                const unsigned mlen,
                const cclean::AdapterIndex &index,
                const bool full_inform = false)
      : AbstractCclean(aligned_output, bed, db, mode, mlen, full_inform),
        index_(index)  {
      if(mode_ == SINGLE_END) checker = new SimpleCleanFunctor;
      if(mode_ == SINGLE_END_Q) checker = new SimpleQualityCleanFunctor;
    }
    virtual ~SimpleClean() { delete checker; }
    virtual Read operator()(const Read &read, bool *ok);

  private:
    const cclean::AdapterIndex &index_;
    AbstractCleanFunctor *checker; // Checks is adapter in read

    // Here goes functors for clean in different modes
    class SimpleCleanFunctor: public AbstractCleanFunctor {
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
    class SimpleQualityCleanFunctor: public AbstractCleanFunctor {
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

#endif /* JOBWRAPPERS_H_ */

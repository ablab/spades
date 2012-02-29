///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Refactored by: Filipe Ribeiro - Nov 2011 - <crdhelp@broadinstitute.org>
//

#ifndef PATHS__FIND_ERRORS_CORE_H
#define PATHS__FIND_ERRORS_CORE_H



#include "kmers/KmerSpectra.h"



class PC_Params
{
public:
  unsigned n_kmers_min;
  unsigned q_high;
  unsigned n_q_high_max_to_lose;
  unsigned q_sum_min_to_win;
  float    ratio_loser;
  set<size_t> i_reads_print;

  PC_Params() : 
    n_kmers_min(6),
    q_high(20),
    n_q_high_max_to_lose(1),
    q_sum_min_to_win(60),
    ratio_loser(0.25),
    i_reads_print()
  {}
};






class EF_Params
{
public:
  int min_Q_to_support;
  int min_n_bases_to_support;
  double max_QSsum_ratio_to_correct;
  double max_QSsum_ratio_to_correct2;
  double max_QSsum_ratio_to_confirm;
  int min_QSsum_to_confirm;
  int min_QSsum_to_win;
  int min_readstack_depth;
  int max_readstack_depth;
  int min_basestack_depth;
  bool skip_under;
  int auto_confirm;
  int print_stack_n;
  int print_stack_n0;
  int print_stack_min_depth;
  int print_stack_max_depth;
  int print_stack_style;
  bool do_palindromes;
  bool do_branches;
  int min_QSsum_to_support;
  int min_n_bases_to_confirm;
  int qual_ceil_radius;
  bool qcr_always;
  set<size_t> i_reads_print;

  EF_Params(int p1, int p2, double p3, double p4, double p5, int p6, int p7, 
            int p8, int p9, int p10, bool p11, int p12, int p13, int p14, 
            int p15, int p16, int p17, bool p18, bool p19, int p20, int p21,
            int p22, bool p23, vec<longlong> plast)
  : min_Q_to_support(p1),
    min_n_bases_to_support(p2),
    max_QSsum_ratio_to_correct(p3),
    max_QSsum_ratio_to_correct2(p4),
    max_QSsum_ratio_to_confirm(p5),
    min_QSsum_to_confirm(p6),
    min_QSsum_to_win(p7),
    min_readstack_depth(p8),
    max_readstack_depth(p9),
    min_basestack_depth(p10),
    skip_under(p11),
    auto_confirm(p12),
    print_stack_n(p13),
    print_stack_n0(p14),
    print_stack_min_depth(p15),
    print_stack_max_depth(p16),
    print_stack_style(p17),
    do_palindromes(p18),
    do_branches(p19),
    min_QSsum_to_support(p20),
    min_n_bases_to_confirm(p21),
    qual_ceil_radius(p22),
    qcr_always(p23)    
  {
    i_reads_print.insert(plast.begin(), plast.end());
  };
};



static 
const EF_Params efp_default(20,               // min_Q_to_support(p1)
                            2,                // min_n_bases_to_support(p2)
                            0.25,             // max_QSsum_ratio_to_correct(p3)
                            0.25,             // max_QSsum_ratio_to_correct2(p4)
                            1.0,              // max_QSsum_ratio_to_confirm(p5)
                            90,               // min_QSsum_to_confirm(p6)
                            40,               // min_QSsum_to_win(p7)
                            5,                // min_readstack_depth(p8)
                            0,                // max_readstack_depth(p9)  0 means infinity
                            1,                // min_basestack_depth(p10)
                            true,             // skip_under(p11)
                            -1,               // auto_confirm(p12)
                            0,                // print_stack_n(p13)
                            0,                // print_stack_n0(p14)
                            5,                // print_stack_min_depth(p15)
                            500,              // print_stack_max_depth(p16)
                            3,                // print_stack_style(p17)
                            false,            // do_palindromes(p18)
                            false,            // do_branches(p19)
                            60,               // min_QSsum_to_support(p20)
                            3,                // min_n_bases_to_confirm(p21)
                            2,                // qual_ceil_radius(p22)
                            true,             // qcr_always(p23)
                            vec<longlong>()); // i_reads_to_print(plast)




class BVLocBase
{
public:
  uint64_t i_read         : 36;   // up to 2^36 = 64G
  uint64_t i_base         : 24;   // up to 2^24 = 16M
  uint64_t rc             :  1;   // for kmers, whether it is RC of FW
  uint64_t palindrome     :  1;   // for kmers, whether it is a palindrome
  uint64_t base           :  2;   // to store a new base
  
  friend
  bool operator < (const BVLocBase & a, const BVLocBase & b)
  { 
    if (a.i_read < b.i_read) return true;
    if (a.i_read > b.i_read) return false;
    return (a.i_base < b.i_base);
  }
};








template<class KMER_t>
class KmerBVLocBase : public KMER_t, public BVLocBase
{
public:
 explicit KmerBVLocBase(const unsigned K = 0) : KMER_t(K) {}

  friend
  bool operator < (const KmerBVLocBase & a, const KmerBVLocBase & b)
  { 
    if (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b)) return true;
    if (static_cast<const KMER_t &>(b) < static_cast<const KMER_t &>(a)) return false;
    return (static_cast<const BVLocBase &>(a) < static_cast<const BVLocBase &>(b));
  }
};















// ---- PreCorrect algorithm
//
// Corrects obvious substitution errors in reads.
// The Precorrect module examines all 25-mers in all reads.  It performs a
// large, external sort in such a way as to bring together for consideration
// all the 25-mers that agree about the initial and final 12-mers, ignoring the
// central base.  For the each pile of 25-mers that agree about the flanking
// bases, and potentially disagree only about that central base, it examines
// the quality scores associated with the base call of the central base.  Piles
// of fewer than 6 are ignored.  For larger piles, the quality scores are
// summed separately for each of the four potential calls of the central base.
// We declare a winner among the four calls as the call having the greatest sum
// of quality scores, but only if the sum is 60 or greater.  If we have a
// winner, then we can also check for losers.  A loser call has no more than
// one call of quality 20 or more, and its sum of quality scores must be less
// than 1/4 that of the winner.  Reads containing loser-calls are candidates for
// having those calls corrected to the winning call.  One further round of
// sorting and screening applies the correction to a read only if it is well
// isolated from other proposed corrections:  no changes are made to either
// call when two proposed corrections are within 12 bases of each other (i.e.,
// when the flanking bases that built a pile are themselves suspect).  Whenever
// a correction is adopted, the associated quality score for that call is set
// to 0.
  

template<class QVV_t>
void pre_correct_parallel(const PC_Params  & pcp, 
                          const unsigned     K,
                          BaseVecVec       * bases_p, 
                          QVV_t            * quals_p, 
                          KmerSpectrum     * kspec_p, 
                          const unsigned     VERBOSITY,
                          const unsigned     NUM_THREADS, 
                          const size_t       mem_mean_ceil = 0);








template<class QVV_t>
void find_errors_parallel(const EF_Params  & efp, 
                          const size_t       K, 
                          BaseVecVec       * bases_p, 
                          QVV_t            * quals_p, 
                          const size_t       NUM_CYCLES,
                          const unsigned     VERBOSITY,
                          const unsigned     NUM_THREADS, 
                          const size_t       mem_mean_ceil = 0);




// ---- Single thread call to find errors
//      Useful when correcting small local read sets in modules other than FindErrors

template<class QVV_t>
void find_errors(const EF_Params  & efp, 
                 const size_t       K, 
                 BaseVecVec       * bases_p, 
                 QVV_t            * quals_p, 
                 const size_t       NUM_CYCLES,
                 const size_t       VERBOSITY);























#endif

////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FirstLookupFinder.h
 * \author tsharpe
 * \date Jun 16, 2009
 *
 * \brief Class that does the "first lookup" operation.
 *
 * Attempts a no-indels alignment of some query sequence by looking up the first kmer
 * of the query (and the first kmer of its reverse complement) in a lookup table, and
 * extending those perfect matches of the first K bases as far as possible.  You can
 * specify filtering criteria for acceptable alignments.
 */
#ifndef LOOKUP_FIRSTLOOKUPFINDERECJ_H_
#define LOOKUP_FIRSTLOOKUPFINDERECJ_H_
#include "feudal/QualNibbleVec.h"
#include "lookup/FirstLookAlign.h"
#include "lookup/FirstLookupFilterECJ.h"
#include "lookup/LookupTab.h"
#include "Basevector.h"
#include "Vec.h"
#include <list>
#include <algorithm>

class FirstLookupFinderECJ
{
private:
  const FirstLookupFilterECJ & _filter;
  const LookupTab & _lookup_tab;
  const BaseVecVec & _contigs;

  const size_t _K;

  // number of reads to process with one trip to the worklist
  static size_t const CHUNK_SIZE = 10000; 


public:
  /// constructor does not copy its args.  they must live as long as the finder.
  FirstLookupFinderECJ(const FirstLookupFilterECJ & filter, 
                       const LookupTab & lookup_tab, 
                       const BaseVecVec & contigs, 
                       const size_t K)
    : _filter(filter), 
      _lookup_tab(lookup_tab), 
      _contigs(contigs), 
      _K(K)
  {}

  // compiler-supplied copy construction and destructor are fine.
  // class is not assignable due to references.

  /// Returns a list of "first lookup" alignments for a specified query.
  /// These are the least-mismatched alignments extended from an initial perfectly matched kmer
  /// for the query or its reverse complement, and subject to the filtering criteria.
  void getAlignments(const BaseVec & query, 
                     const QualNibbleVec & qual, 
                     const uint64_t queryID, 
                     std::list<first_look_align> * result) const;

  /// Version without quality scores
  void getAlignments(const BaseVec & query, 
                     const uint64_t queryID, 
                     std::list<first_look_align> * result) const;

  /// The alignments vector is populated with all alignments for all the queries.
  /// The result will be identical to adding the results of calling getAlignments on each query in turn to
  /// the output vector, but the operation is actually done in parallel, so the results won't be in
  /// exactly the same order as if it had been done that way.
  /// If destructive = true, this destroys queries and quals, to save memory.
  void getAllAlignments(const BaseVecVec & queries, 
                        const VecQualNibbleVec & quals, 
                        vec<first_look_align> * first_aligns, 
                        const size_t n_threads) const;



private:
  void getAlignmentsInternal(const BaseVec & query, 
                     const QualNibbleVec * qualp, 
                     const uint64_t queryID, 
                     std::list<first_look_align> * result) const;

  bool isEligibleF(const Location & location, 
                   const BaseVec &) const
  { 
    size_t nBasesFollowingKmerStart = 
      _contigs[location.getContig()].size() - location.getOffset();
    return (nBasesFollowingKmerStart >= _K && 
            nBasesFollowingKmerStart >= _filter.min_size);
  }

  // For some reason we're rejecting queries that hang off the beginning of a contig,
  // even if _filter.min_size isn't set. [Pretty sure that is so that we always end up
  // with an alignment with offset zero on the query sequence. --bruce]
  bool isEligibleR(const Location & location, 
                   const BaseVec & query) const
  { 
    size_t nBasesPrecedingKmerEnd = location.getOffset() + _lookup_tab.getK();
    size_t nBasesFollowingKmerStart = 
      _contigs[location.getContig()].size() - location.getOffset();
    //return nBasesPrecedingKmerEnd >= query.size() && nBasesFollowingKmerStart >= _filter.min_size;
    //return nBasesPrecedingKmerEnd >= query.size();
    return (nBasesPrecedingKmerEnd >= _K && 
            nBasesPrecedingKmerEnd >= _filter.min_size);
  }

  // count eligible locations for forward alignments
  size_t countEligibleF(const BaseVec & query) const;
  // count eligible locations for reversed alignments
  size_t countEligibleR(const BaseVec & query) const;

  // do forward alignments
  size_t alignF(const BaseVec & query, 
                const uint64_t query_ID) const;
  // do reversed alignments
  size_t alignR(const BaseVec & query, 
                const uint64_t query_ID) const;

  // do forward alignments
  void scoreF(const BaseVec & query, 
              const QualNibbleVec * qual, 
              const uint64_t query_ID, 
              const size_t matchLength, 
              vec<first_look_align> & alignments, 
              vec<unsigned int> & scores) const;
  // do reversed alignments
  void scoreR(const BaseVec & query, 
              const QualNibbleVec * qual, 
              const uint64_t query_ID, 
              const size_t matchLength, 
              vec<first_look_align> & alignments, 
              vec<unsigned int> & scores) const;

  int maxMismatchCount(const size_t queryLen, 
                       const size_t alignmentLen, 
                       const size_t bestMatchCount) const
  { 
    int result = min(static_cast<int>(queryLen * _filter.max_error_rate),
                     static_cast<int>(alignmentLen-bestMatchCount));
    if (_filter.min_match) 
      result = min(result, (int)queryLen - _filter.min_match);
    return result; 
  }

};

// Per-base quality used for scoring when quals not available

#define DEFAULT_QUAL 20



#endif /* LOOKUP_FIRSTLOOKUPFINDERECJ_H_ */

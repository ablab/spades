/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FirstLookupFilter.h
 * \author tsharpe
 * \date Jun 17, 2009
 *
 * \brief Filters alignments produced by the FirstLookup function, or the FirstLookupFinder object.
 */
#ifndef FIRSTLOOKUPFILTERECJ_H_
#define FIRSTLOOKUPFILTERECJ_H_

// FirstLookupFilter
// A struct to hold filtering options for a call to FirstLookup.  These
// filters are all applied early on in the FirstLookup algorithm and serve to
// reduce runtime.  The default behavior is to apply no filters.
struct FirstLookupFilterECJ {

  // Constructor: Sets default values (i.e., no filtering).
  // [Actually, not true...the new scoring parameters do some filtering...
  // either need to clean up documentation or defaults.  --bruce]
  FirstLookupFilterECJ( ) {
    orientation = ALL;
    max_kmer_freq = 0U;
    min_size = 0U;
    max_extend = 0U;
    mismatch_threshhold = 3;
    mismatch_neighborhood = 8;
    mismatch_backoff = 3;
    max_error_rate = 1.f;
    min_match = 0;
    score_max = -1;
    score_delta = 0;
    max_placements = 0;
    vec<int> debug_reads;
  }

  // Filtering options, listed in the order in which they are applied.

  // Only look for alignments where the orientation between query and target is as follows:
  enum { ALL, FW_ONLY, RC_ONLY } orientation;

  // Ignore reads whose leading kmers (using the lookup table's K) appear too
  // often in the lookup table.  This criterion is applied separately to the forward
  // and reverse directions before applying any other criteria.  A zero means that this
  // criterion is not applied.
  unsigned int max_kmer_freq;

  // Only extend from kmers where the potential alignment length is at least this great.
  // This is a criterion on how close to the end (beginning for reverse) of the
  // contig we landed, as well as the query length.
  unsigned int min_size;

  // Maximum number, after filtering by the other criteria above, of potential
  // alignments to try to extend in each orientation.
  // But 0 means that this criterion is not applied.
  // Note that this is a per-chunk number, so the results will differ depending on
  // how you chunk the lookup data.
  // Note: There's a bug in FirstLookup that causes this number to be interpreted as
  // the number to extend, at most.  In other words, rather than skipping expansion
  // if there are too many hits, we always expand the first max_extend hits.
  // FirstLookupFinder faithfully imitates this behavior for regression purposes.
  unsigned int max_extend;

  // These three parameters determine when to quit matching; if it sees mismatch_threshhold
  // incorrect bases in the most recent mismatch_neighborhood, it will then cut off the
  // match at "backoff" bases before the first bad base within the neighborhood.  This is deisgned
  // to identify jumping read junction crossings.
  int mismatch_threshhold;
  int mismatch_neighborhood;
  int mismatch_backoff;

  // If set, require this many bases to align (after cutoff and backoff described by the above
  // parameters) before accepting the alignment.
  int min_match;

  // Max error rate.
  float max_error_rate;

  // Sum of mismatched base qualities much not exceed this.  If negative, don't
  // check.
  int score_max;

  // Score delta: include all aligns with a score within this delta of
  // the best score (i.e., if this value is 0, only those alignment(s)
  // tied for the best score will be returned).  If positive, allows
  // for returning alignments with scores that close to the
  // best.
  unsigned int score_delta;

  // max_placements: if greater than zero, don't bother returning more
  // than this number of placements even if they are close in score.
  unsigned int max_placements;


  // Not really a filtering option, but print debugging information about these reads, if set.
  vec<int> debug_reads;
};

#endif /* FIRSTLOOKUPFILTERECJ_H_ */

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
#ifndef FIRSTLOOKUPFILTER_H_
#define FIRSTLOOKUPFILTER_H_

// FirstLookupFilter
// A struct to hold filtering options for a call to FirstLookup.  These
// filters are all applied early on in the FirstLookup algorithm and serve to
// reduce runtime.  The default behavior is to apply no filters.
struct FirstLookupFilter {

  // Constructor: Sets default values (i.e., no filtering).
  FirstLookupFilter( ) {
    orientation = ALL;
    max_kmer_freq = 0U;
    min_size = 0U;
    max_extend = 0U;
    match_counting_mismatch_threshhold = 4;
    min_match = false;
    max_error_rate = 1.f;
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

  // judging the best match count occurs when we've accumulated this many mismatches
  // setting it to -1 will result in getting simple, highest-identity alignments instead of
  // best identity at the position where a certain mismatch count is achieved.
  int match_counting_mismatch_threshhold;

  // If set, require at least 50% of read to match before accepting an alignment.
  bool min_match;

  // Max error rate.
  float max_error_rate;
};

#endif /* FIRSTLOOKUPFILTER_H_ */

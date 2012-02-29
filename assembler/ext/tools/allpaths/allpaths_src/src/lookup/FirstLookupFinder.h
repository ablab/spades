/////////////////////////////////////////////////////////////////////////////
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
#ifndef LOOKUP_FIRSTLOOKUPFINDER_H_
#define LOOKUP_FIRSTLOOKUPFINDER_H_
#include "lookup/FirstLookupFilter.h"
#include "lookup/FirstLookAlign.h"
#include "lookup/LookupTab.h"
#include "lookup/LookAlign.h"
#include "Basevector.h"
#include "Vec.h"
#include <list>
#include <algorithm>

class FirstLookupFinder
{
public:
    /// constructor does not copy its args.  they must live as long as the finder.
    FirstLookupFinder( FirstLookupFilter const& filter, LookupTab const& lookupTab, vecbvec const& contigs )
    : mFilter(filter), mLookupTab(lookupTab), mContigs(contigs)
    {}

    // compiler-supplied copy construction and destructor are fine.
    // class is not assignable due to references.

    /// Returns a list of "first lookup" alignments for a specified query.
    /// These are the least-mismatched alignments extended from an initial perfectly matched kmer
    /// for the query or its reverse complement, and subject to the filtering criteria.
  void getAlignments( bvec const& query, unsigned int queryID, std::list<first_look_align> & result ) const;

    /// The alignments vector is populated with all alignments for all the queries.
    /// The result will be identical to adding the results of calling getAlignments on each query in turn to
    /// the output vector, but the operation is actually done in parallel, so the results won't be in
    /// exactly the same order as if it had been done that way.
    void getAllAlignments( vecbvec const& queries, vec<look_align>& alignments, unsigned int nThreads = 8 ) const;

private:
    bool isEligibleF( Location const& location, bvec const& ) const
    { unsigned int nBasesFollowingKmerStart = mContigs[location.getContig()].size() - location.getOffset();
      return nBasesFollowingKmerStart >= mFilter.min_size; }

    // For some reason we're rejecting queries that hang off the beginning of a contig,
    // even if mFilter.min_size isn't set.
    // And we're applying the mFilter.min_size criterion to see how close to the end of
    // the contig we are, even though the alignment goes in the other direction.
    // Seems odd, but I believe it's what the old code was doing.
    bool isEligibleR( Location const& location, bvec const& query ) const
    { unsigned int nBasesPrecedingKmerEnd = location.getOffset() + mLookupTab.getK();
      unsigned int nBasesFollowingKmerStart = mContigs[location.getContig()].size() - location.getOffset();
      return nBasesPrecedingKmerEnd >= query.size() && nBasesFollowingKmerStart >= mFilter.min_size; }

    // count eligible locations for forward alignments
    unsigned int countEligibleF( bvec const& query ) const;
    // count eligible locations for reversed alignments
    unsigned int countEligibleR( bvec const& query ) const;

    // do forward alignments
    void alignF( bvec const& query, const int query_ID, std::list<first_look_align>& alignments, unsigned int& bestMatchCount ) const;
    // do reversed alignments
    void alignR( bvec const& query, const int query_ID, std::list<first_look_align>& alignments, unsigned int& bestMatchCount ) const;

    int maxMismatchCount( unsigned int queryLen, unsigned int alignmentLen, unsigned int bestMatchCount ) const
    { int result = min(static_cast<int>(queryLen*mFilter.max_error_rate),static_cast<int>(alignmentLen-bestMatchCount));
      // doesn't actually guarantee 50%, but I believe it's what the old code was doing.
      if ( mFilter.min_match ) result = min(result,static_cast<int>(alignmentLen-(queryLen-4)/2));
      return result; }

    FirstLookupFilter const& mFilter;
    LookupTab const& mLookupTab;
    vecbvec const& mContigs;
    static unsigned int const CHUNK_SIZE = 10000; // number of reads to process with one trip to the worklist
};

#endif /* LOOKUP_FIRSTLOOKUPFINDER_H_ */

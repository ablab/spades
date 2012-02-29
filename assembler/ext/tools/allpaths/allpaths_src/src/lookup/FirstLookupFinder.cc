/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FirstLookupFinder.cc
 * \author tsharpe
 * \date Jun 16, 2009
 *
 * \brief Class that does the "first lookup" operation.
 */
// MakeDepend: library PTHREAD
#include "lookup/FirstLookupFinder.h"
#include "system/Worklist.h"
#include <pthread.h>

#define IMITATE_BUG 1
#define DUMP_BASES 0
#if DUMP_BASES
#include <iostream>
template <class Itr>
void dumpBases( Itr begin, Itr const& end )
{ while ( begin != end )
  { std::cout << Base::val2Char(*begin); ++begin; }
  std::cout << std::endl; }
#endif

namespace
{

class Processor
{
public:
    Processor( FirstLookupFinder const& flf, vecbvec const& queries,
               vec<first_look_align>& alignments, pthread_mutex_t& lock,
               unsigned int chunkSize )
    : mFLF(flf), mQueries(queries), mAlignments(alignments), mLock(lock),
      mChunkSize(chunkSize)
    {}

    Processor( Processor const& that )
    : mFLF(that.mFLF), mQueries(that.mQueries), mAlignments(that.mAlignments),
      mLock(that.mLock), mChunkSize(that.mChunkSize)
    {}

    // copy-assignment prohibited by reference members, compiler-supplied destructor is OK

    void operator()( unsigned int queryID ) const
    {
        unsigned int end = std::min(queryID+mChunkSize,
                                    static_cast<unsigned int>(mQueries.size()));
        std::list<first_look_align> results;
        for ( ; queryID < end; ++queryID )
        {
	  std::list<first_look_align> aligns;
	  mFLF.getAlignments(mQueries[queryID],queryID, aligns );
	  results.insert(results.end(),aligns.begin(),aligns.end());
        }
        if ( results.size() )
        {
            pthread_mutex_lock(&mLock);
            mAlignments.insert(mAlignments.end(),results.begin(),results.end());
            pthread_mutex_unlock(&mLock);
        }
    }

private:
    FirstLookupFinder const& mFLF;
    vecbvec const& mQueries;
    vec<first_look_align>& mAlignments;
    pthread_mutex_t& mLock;
    unsigned int mChunkSize;
};

}

void FirstLookupFinder::getAlignments( bvec const& query, unsigned int queryID,
                                std::list<first_look_align> & result ) const
{
#if DUMP_BASES
std::cout << "Query: " << queryID << std::endl;
#endif
    if ( query.size() >= mLookupTab.getK() && query.size() >= mFilter.min_size )
    {
        unsigned int bestMatchCount = 0;

#if IMITATE_BUG
        alignF(query, queryID, result, bestMatchCount);
        alignR(query, queryID, result, bestMatchCount);
#else
        unsigned int locCountF = countEligibleF(query);
        if ( mFilter.max_extend && locCountF > mFilter.max_extend )
            locCountF = 0;

        unsigned int locCountR = countEligibleR(query);
        if ( mFilter.max_extend && locCountR > mFilter.max_extend )
            locCountR = 0;

	if ( locCountF )
	  alignF(query, queryID, result, bestMatchCount);
	if ( locCountR )
	  alignR(query, queryID, result, bestMatchCount);
#endif
    }
}

void FirstLookupFinder::getAllAlignments( vecbvec const& queries,
                                          vec<look_align>& alignments,
                                          unsigned int nThreads ) const
{
    // Reserve memory for alignments.  We make a heuristic estimate for the
    // upper bound of the alignment frequency per read.
    static const double ALIGNS_PER_READ = 1.25;
    vec<first_look_align> first_aligns;
    first_aligns.reserve(queries.size() * ALIGNS_PER_READ);

    pthread_mutex_t lock;
    pthread_mutex_init(&lock,0);

    if ( true )
    {
        Processor processor(*this, queries, first_aligns, lock, CHUNK_SIZE);
        Worklist<unsigned int,Processor> worklist(processor,nThreads);

        unsigned int nnn = queries.size();
        for ( unsigned int iii = 0; iii < nnn; iii += CHUNK_SIZE )
        {
            if ( worklist.add(iii) > 100000 )
                worklist.waitForEmpty();
        }
        // worklist goes out of scope here, so all the work is known to be done before we kill the mutex
    }

    pthread_mutex_destroy(&lock);
    
    // Convert the first_look_aligns into look_aligns and return them in that format.
    alignments.clear( );
    alignments.resize( first_aligns.size( ) );
    for ( uint ii = 0; ii < first_aligns.size( ); ii++ )
      first_aligns[ii].convert_to_look_align( alignments[ii], queries, mContigs,
                                              mLookupTab.getK( ) );
    
}

unsigned int FirstLookupFinder::countEligibleF( bvec const& query ) const
{
    unsigned int nEligibleLocations = 0;

    if ( mFilter.orientation != FirstLookupFilter::RC_ONLY )
    {
        unsigned int kmer = mLookupTab.getKmer(query.Begin(),
                                               query.Begin(mLookupTab.getK()));
        LocationVec const& locations = mLookupTab.getLocations(kmer);
        if ( !mFilter.max_kmer_freq ||
             static_cast<unsigned int>(locations.size()) <= mFilter.max_kmer_freq )
        {
            if ( !mFilter.min_size )
                nEligibleLocations += locations.size();
            else
            {
                LocationVec::const_iterator end(locations.end());
                LocationVec::const_iterator itr(locations.begin());
                for ( ; itr != end; ++itr )
                {
                    if ( isEligibleF(*itr,query) )
                        nEligibleLocations += 1;
                }
            }
        }
    }

    return nEligibleLocations;
}

unsigned int FirstLookupFinder::countEligibleR( bvec const& query ) const
{
    unsigned int nEligibleLocations = 0;

    if ( mFilter.orientation != FirstLookupFilter::FW_ONLY )
    {
        unsigned int kmer = mLookupTab.getKmer(query.RCBegin(query.size()-mLookupTab.getK()));
        LocationVec const& locations = mLookupTab.getLocations(kmer);
        if ( !mFilter.max_kmer_freq ||
             static_cast<unsigned int>(locations.size()) <= mFilter.max_kmer_freq )
        {
            LocationVec::const_iterator end(locations.end());
            LocationVec::const_iterator itr(locations.begin());
            for ( ; itr != end; ++itr )
            {
                if ( isEligibleR(*itr,query) )
                    nEligibleLocations += 1;
            }
        }
    }

    return nEligibleLocations;
}

void FirstLookupFinder::alignF( bvec const& query, const int query_ID,
                                std::list<first_look_align>& alignments,
                                unsigned int& bestMatchCount ) const
{
    if ( mFilter.orientation != FirstLookupFilter::RC_ONLY )
    {
        unsigned int kmer = mLookupTab.getKmer(query.Begin());
        LocationVec const& locations = mLookupTab.getLocations(kmer);
        if ( !mFilter.max_kmer_freq ||
             static_cast<unsigned int>(locations.size()) <= mFilter.max_kmer_freq )
        {
#if IMITATE_BUG
            LocationVec::const_iterator locsEnd(locations.begin()+min(mFilter.max_extend,(unsigned int)locations.size()));
#else
            LocationVec::const_iterator locsEnd(locations.end());
#endif
            LocationVec::const_iterator locsItr(locations.begin());
            for ( ; locsItr != locsEnd; ++locsItr )
            {
                Location const& location = *locsItr;
                if ( isEligibleF(location,query) )
                {
                    bvec const& contig = mContigs[location.getContig()];
                    unsigned int nBasesFollowingKmerStart = contig.size() - location.getOffset();
                    unsigned int alignmentLen = min(query.size(),nBasesFollowingKmerStart);
#if DUMP_BASES
std::cout << "Location: " << location.getContig() << '.' << location.getOffset() << std::endl;
dumpBases(contig.Begin(location.getOffset()),contig.Begin(location.getOffset()+alignmentLen));
dumpBases(query.Begin(),query.End());
#endif
                    int maxMismatches = maxMismatchCount(query.size(),alignmentLen,bestMatchCount);
                    int nMismatches = 0;
                    unsigned int matchCountAtThreshhold = 0;
                    bvec::const_iterator itr(query.Begin(mLookupTab.getK()));
                    bvec::const_iterator end(query.Begin(alignmentLen));
                    bvec::const_iterator targetItr(contig.Begin(location.getOffset()+mLookupTab.getK()));
                    for ( ; nMismatches <= maxMismatches && itr != end; ++itr,++targetItr )
                    {
                        if ( *itr != *targetItr )
                        {
                            if ( ++nMismatches == mFilter.match_counting_mismatch_threshhold )
                            {
                                matchCountAtThreshhold = alignmentLen - (end - itr) - nMismatches + 1;
                                if ( matchCountAtThreshhold < bestMatchCount )
                                    break;
                            }
                        }
                    }
                    if ( nMismatches <= maxMismatches )
                    {
                        if ( nMismatches < mFilter.match_counting_mismatch_threshhold )
                            matchCountAtThreshhold = query.size() - nMismatches;

                        if ( matchCountAtThreshhold >= bestMatchCount )
                        {
                            if ( matchCountAtThreshhold > bestMatchCount )
                            {
                                bestMatchCount = matchCountAtThreshhold;
                                alignments.clear();
                            }

			    // Found an alignment!  Record a first_look_align.
			    alignments.insert( alignments.end( ), 1, first_look_align( ) );
			    alignments.back( ).target_loc = location;
			    alignments.back( ).query_ID = query_ID;
			    alignments.back( ).n_mismatches = nMismatches;
                        }
                    }
                }
            }
        }
    }
}

void FirstLookupFinder::alignR( bvec const& query, const int query_ID,
                                std::list<first_look_align>& alignments,
                                unsigned int& bestMatchCount ) const
{
    if ( mFilter.orientation != FirstLookupFilter::FW_ONLY )
    {
        unsigned int kStart = query.size() - mLookupTab.getK();
        unsigned int kmer = mLookupTab.getKmer(query.RCBegin(kStart));
        LocationVec const& locations = mLookupTab.getLocations(kmer);
        if ( !mFilter.max_kmer_freq ||
             static_cast<unsigned int>(locations.size()) <= mFilter.max_kmer_freq )
        {
#if IMITATE_BUG
            LocationVec::const_iterator locsEnd(locations.begin()+min(mFilter.max_extend,(unsigned int)locations.size()));
#else
            LocationVec::const_iterator locsEnd(locations.end());
#endif
            LocationVec::const_iterator locsItr(locations.begin());
            for ( ; locsItr != locsEnd; ++locsItr )
            {
                Location const& location = *locsItr;
                if ( isEligibleR(location,query) )
                {
                    bvec const& contig = mContigs[location.getContig()];
                    unsigned int nBasesPrecedingKmerEnd = location.getOffset() + mLookupTab.getK();
                    unsigned int alignmentLen = min(query.size(),nBasesPrecedingKmerEnd);
                    unsigned int offset = nBasesPrecedingKmerEnd - alignmentLen;
#if DUMP_BASES
std::cout << "Location: " << location.getContig() << '.' << offset << std::endl;
dumpBases(contig.Begin(offset),contig.Begin(offset+alignmentLen));
dumpBases(query.RCBegin(query.size()-alignmentLen),query.RCEnd());
#endif
                    int maxMismatches = maxMismatchCount(query.size(),alignmentLen,bestMatchCount);
                    int nMismatches = 0;
                    unsigned int matchCountAtThreshhold = 0;
                    bvec::const_iterator itr(query.Begin(mLookupTab.getK()));
                    bvec::const_iterator end(query.Begin(alignmentLen));
                    bvec::const_rc_iterator targetItr(contig.RCBegin(contig.size()-location.getOffset()));
                    for ( ; nMismatches <= maxMismatches && itr != end; ++itr,++targetItr )
                    {
                        if ( *itr != *targetItr )
                        {
                            if ( ++nMismatches == mFilter.match_counting_mismatch_threshhold )
                            {
                                matchCountAtThreshhold = alignmentLen - (end - itr) - nMismatches + 1;
                                if ( matchCountAtThreshhold < bestMatchCount )
                                    break;
                            }
                        }
                    }
                    if ( nMismatches <= maxMismatches )
                    {
                        if ( nMismatches < mFilter.match_counting_mismatch_threshhold )
                            matchCountAtThreshhold = query.size() - nMismatches;

                        if ( matchCountAtThreshhold >= bestMatchCount )
                        {
                            if ( matchCountAtThreshhold > bestMatchCount )
                            {
                                bestMatchCount = matchCountAtThreshhold;
                                alignments.clear();
                            }

			    // Found an alignment!  Record a first_look_align.
			    alignments.insert( alignments.end( ), 1, first_look_align( ) );
			    alignments.back( ).target_loc = location;
			    alignments.back( ).query_ID = query_ID;
			    alignments.back( ).n_mismatches = -nMismatches - 1; // flip sign of n_mismatches to indicate RC orientation
                        }
                    }
                }
            }
        }
    }
}

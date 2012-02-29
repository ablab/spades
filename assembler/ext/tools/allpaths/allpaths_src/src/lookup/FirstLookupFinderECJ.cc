/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file FirstLookupFinderECJ.cc
 * \author tsharpe
 * \date Jun 16, 2009
 *
 * \brief Class that does the "first lookup" operation.
 */
// MakeDepend: library PTHREAD
#include "VecUtilities.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "system/Worklist.h"
#include <pthread.h>

#define IMITATE_BUG 1


#define DUMP_BASES 0
#if DUMP_BASES
#include <iostream>
template <class Itr>
void dumpBases(Itr begin, Itr const& end)
{ 
  while (begin != end) { 
    std::cout << Base::val2Char(*begin); ++begin; 
  }
  std::cout << std::endl; 
}
#endif



#define is_in(x,v) (find(v.begin(), v.end(), x) != v.end())



namespace
{
  class AlignProcessor
  {
  private:
    const FirstLookupFinderECJ & _flf;
    const BaseVecVec           & _queries;
    const VecQualNibbleVec     & _quals;

    vec<first_look_align>      * _first_aligns;

    pthread_mutex_t            & _lock;
    const size_t                 _chunk_size;

  public:
    AlignProcessor(const FirstLookupFinderECJ & flf, 
                   const BaseVecVec & queries, 
                   const VecQualNibbleVec & quals,
                   vec<first_look_align> * first_aligns, 
                   pthread_mutex_t& lock,
                   const size_t chunkSize)
      : _flf(flf), 
        _queries(queries), 
        _quals(quals), 
        _first_aligns(first_aligns), 
        _lock(lock),
	_chunk_size(chunkSize)
    {}

    AlignProcessor(const AlignProcessor & that)
      : _flf(that._flf), 
        _queries(that._queries), 
        _quals(that._quals), 
        _first_aligns(that._first_aligns),
	_lock(that._lock), 
        _chunk_size(that._chunk_size)
    {}

    // copy-assignment prohibited by reference members
    // compiler-supplied destructor is OK

    void operator()(size_t query_ID) const
    {
      size_t end = std::min(query_ID + _chunk_size, _queries.size());

      std::list<first_look_align> results;

      for (; query_ID < end; ++query_ID) {
        std::list<first_look_align> aligns;
        _flf.getAlignments(_queries[query_ID], _quals[query_ID], query_ID, & aligns);
        results.insert(results.end(), aligns.begin(), aligns.end());
      }

      if (results.size()) {
        pthread_mutex_lock(&_lock);
        _first_aligns->insert(_first_aligns->end(), results.begin(), results.end());
        pthread_mutex_unlock(&_lock);
      }
    }
   
  };
}



void FirstLookupFinderECJ::getAlignmentsInternal(const BaseVec & query, 
						 const QualNibbleVec * qual,
						 const uint64_t query_ID,
						 std::list<first_look_align> * result) const
{
#if DUMP_BASES
  std::cout << "Query: " << query_ID << std::endl;
#endif
  vec<first_look_align> aligns; 
  if (query.size() >= _lookup_tab.getK() && 
      query.size() >= _filter.min_size) {
    
    const bool debug = is_in(query_ID, _filter.debug_reads);



    // Align[FR] go through all the alignments and determine the
    // maximum matching length (searching for potential jump read
    // junction crossing).  Score[FR] then score each candidate
    // alignment up through that length using quality scores. --bruce

    size_t bestMatchCount = 0;
#if IMITATE_BUG

    bestMatchCount = max(alignF(query, query_ID), bestMatchCount);
    if (debug) PRINT2(query_ID, bestMatchCount);
    bestMatchCount = max(alignR(query, query_ID), bestMatchCount);
    if (debug) PRINT2(query_ID, bestMatchCount);

    if (_filter.min_match > 0 && 
        bestMatchCount < static_cast<unsigned>(_filter.min_match)) 
      return;

    vec<unsigned int> scores;

    scoreF(query, qual, query_ID, bestMatchCount, aligns, scores);
    if (debug) PRINT2(query_ID, scores.size());
    scoreR(query, qual, query_ID, bestMatchCount, aligns, scores);
    if (debug) PRINT2(query_ID, scores.size());

#else

    size_t locCountF = countEligibleF(query);
    if (_filter.max_extend && locCountF > _filter.max_extend)
      locCountF = 0;

    size_t locCountR = countEligibleR(query);
    if (_filter.max_extend && locCountR > _filter.max_extend)
      locCountR = 0;


    if (locCountF)
      bestMatchCount = max(alignF(query, query_ID), bestMatchCount);
    if (locCountR)
      bestMatchCount = max(alignR(query, query_ID), bestMatchCount);

    if (_filter.min_match > 0 && 
        bestMatchCount < _filter.min_match) 
      return;

    if (locCountF)
      scoreF(query, qual, query_ID, bestMatchCount, aligns, scores);
    if (locCountR)
      scoreR(query, qual, query_ID, bestMatchCount, aligns, scores);

#endif

    if (scores.size() > 0) {

      SortSync(scores, aligns);

      size_t n_to_return = scores.size();
      if (_filter.max_placements > 0 && 
          _filter.max_placements < n_to_return)
        n_to_return = _filter.max_placements;

      size_t bestscore = scores[0];
      if (debug) PRINT5(query_ID, _filter.score_max, scores.size(), bestscore, n_to_return);
      if (_filter.score_max < 0 || bestscore <= static_cast<unsigned>(_filter.score_max)) {
        // if there is a perfect score, don't return imperfect alternatives
        size_t score_limit = (bestscore==0) ? 0 : (bestscore + _filter.score_delta);
        for (size_t i = 0; i < n_to_return && scores[i] <= score_limit; ++i) {

          result->insert(result->end(), 1, aligns[i]);

          if (debug) {
            PRINT3(scores[i], 
                   aligns[i].target_loc.getContig(), 
                   aligns[i].target_loc.getOffset());
          }
        }

      }

    }

  }
}


void FirstLookupFinderECJ::getAlignments(const BaseVec & query, 
						 const QualNibbleVec & qual,
						 const uint64_t query_ID,
						 std::list<first_look_align> * result) const
{
  getAlignmentsInternal(query, &qual, query_ID, result);
}

void FirstLookupFinderECJ::getAlignments(const BaseVec & query, 
                                         const uint64_t query_ID,
                                         std::list<first_look_align> * result) const
{
  getAlignmentsInternal(query, NULL, query_ID, result);
}


void FirstLookupFinderECJ::getAllAlignments(const BaseVecVec & queries,
                                            const VecQualNibbleVec & quals,
                                            vec<first_look_align> * first_aligns,
                                            const size_t n_threads) const
{
  // Reserve memory for alignments.  We make a heuristic estimate for the
  // upper bound of the alignment frequency per read.
  static const double ALIGNS_PER_READ = 1.25;
  first_aligns->reserve(queries.size() * ALIGNS_PER_READ);

  pthread_mutex_t lock;
  pthread_mutex_init(&lock, 0);

  if (true) {
    AlignProcessor processor(*this, queries, quals, 
                             first_aligns, 
                             lock, CHUNK_SIZE);
    Worklist<size_t, AlignProcessor> worklist(processor, n_threads - 1);
    
    size_t nnn = queries.size();
    for (size_t iii = 0; iii < nnn; iii += CHUNK_SIZE) {
      if (worklist.add(iii) > 100000)
        worklist.waitForEmpty();
    }
    // worklist goes out of scope here, 
    // so all the work is known to be done before we kill the mutex
  }
 
  pthread_mutex_destroy(&lock);
}








size_t FirstLookupFinderECJ::countEligibleF(const BaseVec& query) const
{
  size_t nEligibleLocations = 0;

  if (_filter.orientation != FirstLookupFilterECJ::RC_ONLY) {
    size_t kmer = _lookup_tab.getKmer(query.Begin(),
                                      query.Begin(_lookup_tab.getK()));
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
      if (!_filter.min_size)
        nEligibleLocations += locations.size();
      else {
        LocationVec::const_iterator end(locations.end());
        LocationVec::const_iterator itr(locations.begin());
        for (; itr != end; ++itr) {
          if (isEligibleF(*itr,query))
            nEligibleLocations += 1;
        }
      }
    }
  }

  return nEligibleLocations;
}

size_t FirstLookupFinderECJ::countEligibleR(const BaseVec& query) const
{
  size_t nEligibleLocations = 0;

  if (_filter.orientation != FirstLookupFilterECJ::FW_ONLY) {
    size_t kmer = _lookup_tab.getKmer(query.RCBegin(query.size()-_lookup_tab.getK()));
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
      LocationVec::const_iterator end(locations.end());
      LocationVec::const_iterator itr(locations.begin());
      for (; itr != end; ++itr) {
        if (isEligibleR(*itr,query))
          nEligibleLocations += 1;
      }
    }
  }

  return nEligibleLocations;
}

size_t FirstLookupFinderECJ::alignF(const BaseVec & query,
                                    const size_t query_ID) const
{
  size_t matchLength = 0;
  bool debug = is_in(query_ID, _filter.debug_reads);
  if (_filter.orientation != FirstLookupFilterECJ::RC_ONLY) {

    unsigned int kmer = _lookup_tab.getKmer(query.Begin());
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    if (debug) PRINT2(query_ID, locations.size());
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
#if IMITATE_BUG
      LocationVec::const_iterator locsEnd(locations.begin()+min(_filter.max_extend,(unsigned int)locations.size()));
#else
      LocationVec::const_iterator locsEnd(locations.end());
#endif
      LocationVec::const_iterator locsItr(locations.begin());
      for (; locsItr != locsEnd; ++locsItr) {
        if (matchLength == query.size()) break;
        Location const& location = *locsItr;
        if (isEligibleF(location,query)) {
          const BaseVec& contig = _contigs[location.getContig()];
          unsigned int nBasesFollowingKmerStart = contig.size() - location.getOffset();
          unsigned int alignmentLen = min(query.size(),nBasesFollowingKmerStart);
#if DUMP_BASES
          std::cout << "Location: " << location.getContig() << '.' << location.getOffset() << std::endl;
          dumpBases(contig.Begin(location.getOffset()),contig.Begin(location.getOffset()+alignmentLen));
          dumpBases(query.Begin(),query.End());
#endif
          bvec::const_iterator itr(query.Begin(_lookup_tab.getK()));
          bvec::const_iterator end(query.Begin(alignmentLen));
          bvec::const_iterator targetItr(contig.Begin(location.getOffset()+_lookup_tab.getK()));
          // Circular memory of recent (mis)matches and index into it.  If not specified,
          // use the whole length.
          int nhood_size = _filter.mismatch_neighborhood;
          if (nhood_size == 0) nhood_size = alignmentLen;
          vec<int> nhood_mismatch(nhood_size, 0);
          int nhood_index = 0, nhood_mismatches = 0;
          int nMismatches = 0;
          unsigned int matchCountAtThreshhold = alignmentLen;
          for (; itr != end; ++itr,++targetItr) {
            // "was": mismatch nood_size bases ago, "is": mismatch this base
            int was = nhood_mismatch[nhood_index];
            int is = (*itr != *targetItr) ? 1 : 0;
            nhood_mismatch[nhood_index] = is;
            // advance circular index
            if (++nhood_index == nhood_size) nhood_index = 0;
            // keep running tally of mismatches within recent neighborhood
            nhood_mismatches += is - was;
            // keep running tally of mismatches we've passed over
            nMismatches += was;
            // if we hit our mismatchlimit in nhood
            if (nhood_mismatches == _filter.mismatch_threshhold) {
              matchCountAtThreshhold = alignmentLen - (end - itr) + 1;
              // go back before first mismatch within neighborhood
              // nhood_index is pointing at least recent slot
              for (int n = 0; n < nhood_size; ++n) {
                int ni = nhood_index + n;
                if (ni >= nhood_size) ni -= nhood_size;
                if (nhood_mismatch[ni]) {
                  matchCountAtThreshhold -= nhood_size - n;
                  break;
                }
              }
              matchCountAtThreshhold -= _filter.mismatch_backoff;
              break;
            }
          }
          // Matches prior to entering the recent neighborhood
          nMismatches = query.size() - matchCountAtThreshhold;
          ForceAssertGe(nMismatches, 0);

          if (matchCountAtThreshhold > matchLength)
            matchLength = matchCountAtThreshhold;
        }
      }
    }
  }
  
  return matchLength;
}

size_t FirstLookupFinderECJ::alignR(const BaseVec & query,
                                    const size_t query_ID) const
{
  size_t matchLength = 0;
  bool debug = is_in(query_ID, _filter.debug_reads);
  if (_filter.orientation != FirstLookupFilterECJ::FW_ONLY) {
    unsigned int kStart = query.size() - _lookup_tab.getK();
    unsigned int kmer = _lookup_tab.getKmer(query.RCBegin(kStart));
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    if (debug) PRINT2(query_ID, locations.size());
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
#if IMITATE_BUG
      LocationVec::const_iterator locsEnd(locations.begin()+min(_filter.max_extend,(unsigned int)locations.size()));
#else
      LocationVec::const_iterator locsEnd(locations.end());
#endif
      LocationVec::const_iterator locsItr(locations.begin());
      for (; locsItr != locsEnd; ++locsItr) {
        if (matchLength == query.size()) break;
        Location const& location = *locsItr;
        if (isEligibleR(location,query)) {
          const BaseVec& contig = _contigs[location.getContig()];
          unsigned int nBasesPrecedingKmerEnd = location.getOffset() + _lookup_tab.getK();
          unsigned int alignmentLen = min(query.size(),nBasesPrecedingKmerEnd);
          unsigned int offset = nBasesPrecedingKmerEnd - alignmentLen;
#if DUMP_BASES
          std::cout << "Location: " << location.getContig() << '.' << offset << std::endl;
          dumpBases(contig.Begin(offset),contig.Begin(offset+alignmentLen));
          dumpBases(query.RCBegin(query.size()-alignmentLen),query.RCEnd());
#endif
          bvec::const_iterator itr(query.Begin(_lookup_tab.getK()));
          bvec::const_iterator end(query.Begin(alignmentLen));
          bvec::const_rc_iterator targetItr(contig.RCBegin(contig.size()-location.getOffset()));
          // Circular memory of recent (mis)matches and index into it.  If not specified,
          // use the whole length.
          int nhood_size = _filter.mismatch_neighborhood;
          if (nhood_size == 0) nhood_size = alignmentLen;
          vec<int> nhood_mismatch(nhood_size, 0);
          int nhood_index = 0, nhood_mismatches = 0;
          int nMismatches = 0;
          unsigned int matchCountAtThreshhold = alignmentLen;
          for (; itr != end; ++itr,++targetItr) {
            // "was": mismatch nood_size bases ago, "is": mismatch this base
            int was = nhood_mismatch[nhood_index];
            int is = (*itr != *targetItr) ? 1 : 0;
		      
            nhood_mismatch[nhood_index] = is;
            // advance circular index
            if (++nhood_index == nhood_size) nhood_index = 0;
            // keep running tally of mismatches within recent neighborhood
            nhood_mismatches += is - was;
            // keep running tally of mismatches we've passed over
            nMismatches += was;
            // if we hit our mismatchlimit in nhood
            if (nhood_mismatches == _filter.mismatch_threshhold) {
              matchCountAtThreshhold = alignmentLen - (end - itr) + 1;
              // go back before first mismatch within neighborhood
              // nhood_index is pointing at least recent slot
              for (int n = 0; n < nhood_size; ++n) {
                int ni = nhood_index + n;
                if (ni >= nhood_size) ni -= nhood_size;
                if (nhood_mismatch[ni]) {
                  matchCountAtThreshhold -= nhood_size - n;
                  break;
                }
              }
              matchCountAtThreshhold -= _filter.mismatch_backoff;
              break;
            }
          }
          nMismatches = query.size() - matchCountAtThreshhold;
          ForceAssertGe(nMismatches, 0);

          if (matchCountAtThreshhold > matchLength)
            matchLength = matchCountAtThreshhold;
        }
      }
    }
  }
  return matchLength;
}


void FirstLookupFinderECJ::scoreF(const BaseVec & query, 
                                  const QualNibbleVec * qual,
                                  const size_t query_ID,
                                  const size_t matchLength,
                                  vec<first_look_align> & alignments,
                                  vec<unsigned int> & scores) const
{
    
  if (_filter.orientation != FirstLookupFilterECJ::RC_ONLY &&
      matchLength > _lookup_tab.getK()) {
    unsigned int kmer = _lookup_tab.getKmer(query.Begin());
    unsigned int perfects = 0;
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    bool debug = is_in(query_ID, _filter.debug_reads);
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
#if IMITATE_BUG
      LocationVec::const_iterator locsEnd(locations.begin()+min(_filter.max_extend,(unsigned int)locations.size()));
#else
      LocationVec::const_iterator locsEnd(locations.end());
#endif
      LocationVec::const_iterator locsItr(locations.begin());
      for (; locsItr != locsEnd; ++locsItr) {
        Location const& location = *locsItr;
        if (debug) {
          PRINT4(query_ID, location.getContig(), _contigs[location.getContig()].size(), location.getOffset());
        }
        if (isEligibleF(location,query)) {
          const BaseVec& contig = _contigs[location.getContig()];
          unsigned int nBasesFollowingKmerStart = contig.size() - location.getOffset();
          unsigned int alignmentLen = min(query.size(),nBasesFollowingKmerStart);
          // If we don't have enough continuity to match the best length, punt..
          if (alignmentLen < matchLength) continue;
          if (debug) PRINT4(query_ID, nBasesFollowingKmerStart, alignmentLen, matchLength);
#if DUMP_BASES
          std::cout << "Location: " << location.getContig() << '.' << location.getOffset() << std::endl;
          dumpBases(contig.Begin(location.getOffset()),contig.Begin(location.getOffset()+alignmentLen));
          dumpBases(query.Begin(),query.End());
#endif
          int nMismatches = 0;
          unsigned int matchCountAtThreshhold = 0;
          //bvec::const_iterator itr(query.Begin(_lookup_tab.getK()));
          //QualNibbleVec::const_iterator qitr((*qual).begin(_lookup_tab.getK()));
          //bvec::const_iterator end(query.Begin(matchLength));
          bvec::const_iterator targetItr(contig.Begin(location.getOffset()+_lookup_tab.getK()));
          int score = 0;
          for (u_int i = _lookup_tab.getK(); i < matchLength; ++i,++targetItr) {
            if (query[i] != *targetItr) {
              score += (qual==NULL) ? DEFAULT_QUAL : (*qual)[i];
              //score += ((qual==NULL) ? DEFAULT_QUAL : (*qual)[i]) / 3;
              ++nMismatches;
              if (debug) cout << "*";
            } else if (debug) cout << ".";
          }
          if (debug) {
            cout << endl;
            PRINT4(query_ID, score, nMismatches, matchLength);
          }
          matchCountAtThreshhold = matchLength - nMismatches;
          nMismatches = query.size() - matchCountAtThreshhold;

          ForceAssertGe(nMismatches, 0);
          ForceAssertGe(score, 0);

          //  Record a first_look_align and its score.
          first_look_align fla;
          fla.target_loc = location;
          fla.query_ID = query_ID;
          fla.n_mismatches = nMismatches;
          alignments.push_back(fla);
          scores.push_back(score);
          if (score == 0 && _filter.max_placements && ++perfects >= _filter.max_placements)
            break;
        }
      }
    }
  }
}

void FirstLookupFinderECJ::scoreR(const BaseVec & query, 
                                  const QualNibbleVec * qual,
                                  const size_t query_ID,
                                  const size_t matchLength,
                                  vec<first_look_align>& alignments,
                                  vec<unsigned int>& scores) const
{
  if (_filter.orientation != FirstLookupFilterECJ::FW_ONLY &&
      matchLength > _lookup_tab.getK()) {
    unsigned int kStart = query.size() - _lookup_tab.getK();
    unsigned int kmer = _lookup_tab.getKmer(query.RCBegin(kStart));
    unsigned int perfects = 0;
    LocationVec const& locations = _lookup_tab.getLocations(kmer);
    bool debug = is_in(query_ID, _filter.debug_reads);
    if (!_filter.max_kmer_freq ||
        static_cast<unsigned int>(locations.size()) <= _filter.max_kmer_freq) {
#if IMITATE_BUG
      LocationVec::const_iterator locsEnd(locations.begin()+min(_filter.max_extend,(unsigned int)locations.size()));
#else
      LocationVec::const_iterator locsEnd(locations.end());
#endif
      LocationVec::const_iterator locsItr(locations.begin());
      for (; locsItr != locsEnd; ++locsItr) {
        Location const& location = *locsItr;
        if (debug) {
          PRINT4(query_ID, location.getContig(), _contigs[location.getContig()].size(), location.getOffset());
        }
        if (isEligibleR(location,query)) {
          const BaseVec& contig = _contigs[location.getContig()];
          unsigned int nBasesPrecedingKmerEnd = location.getOffset() + _lookup_tab.getK();
          unsigned int alignmentLen = min(query.size(),nBasesPrecedingKmerEnd);
          // If we don't have enough continuity to match the best length, punt..
          if (alignmentLen < matchLength) continue;
          unsigned int offset = nBasesPrecedingKmerEnd - alignmentLen;
          if (debug) PRINT3(nBasesPrecedingKmerEnd, alignmentLen, matchLength);
#if DUMP_BASES
          std::cout << "Location: " << location.getContig() << '.' << offset << std::endl;
          dumpBases(contig.Begin(offset),contig.Begin(offset+alignmentLen));
          dumpBases(query.RCBegin(query.size()-alignmentLen),query.RCEnd());
#endif
          //int maxMismatches = maxMismatchCount(query.size(),alignmentLen,matchLength);
          int nMismatches = 0;
          unsigned int matchCountAtThreshhold = 0;
          //bvec::const_iterator itr(query.Begin(_lookup_tab.getK()));
          //QualNibbleVec::const_iterator qitr((*qual).begin(_lookup_tab.getK()));
          //bvec::const_iterator end(query.Begin(matchLength));
          bvec::const_rc_iterator targetItr(contig.RCBegin(contig.size()-location.getOffset()));
          int score = 0;
          for (u_int i = _lookup_tab.getK(); i < matchLength; ++i,++targetItr) {
            if (query[i] != *targetItr) {
              score += (qual==NULL) ? DEFAULT_QUAL : (*qual)[i];
              //score += ((qual==NULL) ? DEFAULT_QUAL : (*qual)[i]) / 3;
              ++nMismatches;
              if (debug) cout << "*";
            } else if (debug) cout << ".";
          }
          if (debug) {
            cout << endl;
            PRINT4(query_ID, score, nMismatches, matchLength);
          }
          matchCountAtThreshhold = matchLength - nMismatches;
          nMismatches = query.size() - matchCountAtThreshhold;

          ForceAssertGe(nMismatches, 0);
          ForceAssertGe(score, 0);

          //  Record a first_look_align and its score.
          first_look_align fla;
          fla.target_loc = location;
          fla.query_ID = query_ID;
          fla.n_mismatches = -nMismatches - 1; // flip sign of n_mismatches to indicate RC orientation
          alignments.push_back(fla);
          scores.push_back(score);
          if (score == 0 && _filter.max_placements && ++perfects >= _filter.max_placements)
            break;
        }
      }
    }
  }
}



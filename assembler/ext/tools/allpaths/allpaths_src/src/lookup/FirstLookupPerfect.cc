/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
// MakeDepend: library PTHREAD
#include "VecUtilities.h"
#include "lookup/FirstLookupPerfect.h"
#include "system/Worklist.h"
#include <pthread.h>

namespace 
{
  class Processor 
  {
  private:
    const FirstLookupPerfect & _flp;
    const BaseVecVec         & _queries;

    vec<first_look_align>    * _alignments;

    pthread_mutex_t          & _lock;
    const size_t               _chunk_size;
    
  public:

    Processor(const FirstLookupPerfect & flp,
              const BaseVecVec & queries,
              vec<first_look_align> * alignments,
              pthread_mutex_t & lock,
              const size_t chunk_size) 
      : _flp(flp),
        _queries(queries),
        _alignments(alignments),
        _lock(lock),
        _chunk_size(chunk_size)
    {}
    
    Processor(const Processor &that) :
      _flp(that._flp),
      _queries(that._queries),
      _alignments(that._alignments),
      _lock(that._lock),
      _chunk_size(that._chunk_size)
    {}
    
    // copy-assignment prohibited by reference members,
    // compiler-supplied destructor is OK
    
    void operator() (size_t query_ID) const 
    {
      size_t end = std::min(query_ID + _chunk_size, _queries.size());

      vec<first_look_align> results;

      for (; query_ID < end; ++query_ID) {
	vec<first_look_align> aligns;
	_flp.getAlignments(_queries[query_ID], query_ID, & aligns);
	copy(aligns.begin(), aligns.end(), back_inserter(results));
      }
      if (results.size()) {
	pthread_mutex_lock(&_lock);
	_alignments->insert(_alignments->end(), results.begin(), results.end());
	pthread_mutex_unlock(&_lock);
      }
    }
    
  };
  
}







void FirstLookupPerfect::getAlignments(const BaseVec & query,
                                       const uint64_t query_ID,
                                       vec<first_look_align> * result) const
{
  vec<first_look_align> aligns; 
  if (query.size() < _lookup_tab.getK()) return;
  
  // Select kmer with lowest frequency.
  const size_t K = _lookup_tab.getK();
  
  uint low_freq = 0;
  uint low_pos = 0;
  uint pos = 0;
  while(pos + K < query.size()) {
    uint kmer = _lookup_tab.getKmer(query.Begin(pos));  
    uint freq = _lookup_tab.getFreq(kmer);

    if (pos == 0 || freq < low_freq) {
      low_pos = pos;
      low_freq = freq;
    }
    pos += K;
  }
  const uint kmer = _lookup_tab.getKmer(query.Begin(low_pos));
  const LocationVec &locations = _lookup_tab.getLocations(kmer);
  uint nlocs = locations.size();
  
  if (nlocs > MAX_FREQ) return;
  
  // Loop over all locations (for the kmer with lowest frequency).
  LocationVec::const_iterator locsEnd(locations.end());
  LocationVec::const_iterator locsItr(locations.begin());
  for (; locsItr != locsEnd; ++locsItr) {
    const Location &location = *locsItr;
    const BaseVec &contig = _contigs[location.getContig()];

    // The read is not completely embedded in the contig.
    if (low_pos > location.getOffset()) continue;
    uint true_offset = location.getOffset() - low_pos;
    if (true_offset + query.size() > contig.size()) continue;
    
    // Match bases.
    BaseVec::const_iterator itr(query.Begin());
    BaseVec::const_iterator end(query.Begin(query.size()));
    BaseVec::const_iterator targetItr(contig.Begin(true_offset));
    bool is_bad = false;
    for (; itr < end; ++itr,++targetItr) {
      if (*itr != *targetItr) {
	is_bad = true;
	break;
      }
    }
    if (is_bad) continue;
    
    first_look_align fla;
    fla.target_loc = location;
    fla.query_ID = query_ID;
    fla.n_mismatches = 0;
    aligns.push_back(fla);

    if (aligns.size() >= MAX_HITS) break;
  }
  
  // Append aligns to results.
  copy(aligns.begin(), aligns.end(), back_inserter(*result));
}

void FirstLookupPerfect::getAllAlignments(const BaseVecVec & queries,
                                          vec<first_look_align> * first_aligns,
                                          const size_t n_threads) const
{
  first_aligns->clear();




#ifdef SKIP
  { // SANTEMP - ABSOLUTELY DO NOT COMMIT THIS

    first_aligns->reserve(queries.size() * 1.5);
    const size_t K = _lookup_tab.getK();
    
    size_t dotter = 1000;
    for (size_t ii = 0; ii < queries.size(); ii++) {
      DotMod(cout, ii, dotter);
      
      vec<first_look_align> fla;
      this->getAlignments(queries[ii], ii, fla);
      
      for (size_t jj = 0; jj < fla.size(); jj++) {
	first_aligns->push_back(fla[jj]);
      }

    }
    
    return;
  } // SANTEMP - ABSOLUTELY DO NOT COMMIT THIS
#endif
  
  
  
  // HEURISTICS - Reserve memory for alignments (estimate upper bound).
  static const double ALIGNS_PER_READ = 1.25;

  // Reserve memory.
  const size_t K = _lookup_tab.getK();
  first_aligns->reserve(queries.size() * ALIGNS_PER_READ);
  
  pthread_mutex_t lock;
  pthread_mutex_init(&lock,0);
  
  if (true) {
    Processor processor(*this, queries, 
                        first_aligns, 
                        lock, CHUNK_SIZE);
    Worklist<size_t, Processor> worklist(processor, n_threads - 1);
    
    size_t nnn = queries.size();
    for (size_t iii = 0; iii < nnn; iii += CHUNK_SIZE) {
      if (worklist.add(iii) > 100000)
	worklist.waitForEmpty();
    }

    // The worklist goes out of scope here, so all the work is known
    // to be done before we kill the mutex.
  }
  
  pthread_mutex_destroy(&lock);
}


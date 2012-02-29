#ifndef LOOKUP_TOOLS_H
#define LOOKUP_TOOLS_H
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Basevector.h"
#include "lookup/KmerIndex.h"
#include "lookup/LookupTable.h"
#include "lookup/PerfectCount.h"

#ifdef __GNUC__
#include <ext/slist>
using __gnu_cxx::slist;
#endif

/// \file LookupTools.h - declarations and template definitions of global methods
/// used to work with lookup tables, queries and hits.

/// Looks ups queries specified by the iterator interval
/// [q_start, q_end) in the lookup table \c lookup_table
/// and records all hits. NOTE: this method retrieves hits
/// <em>from the currently loaded chunk</em> of the lookup table only!
///
/// @param lookup reference genome represented as lookup table
/// @param q_start start of the range of queries
/// @param q_end   end of the range of queries
/// @param[out] hits all hits (occurences on the genome) found for 
/// all queries in the range [q_start, q_end) will be \em appended
/// to this vector (previous content is \em not destroyed).
/// @return total number of hits found.
template <class QueryIterator>
void QueriesToHits(const lookup_table & lookup, 
		   QueryIterator q_start, 
		   QueryIterator q_end, 
		   vec<RawHit> & hits) 
{
    unsigned int nhits = 0;
    for ( ; q_start != q_end ; ++q_start ) {
        // get the range in the lookup table index spanning
        // all the offsets for the current Kmer (query):
 	lookup_table::locs_iterator iloc_begin = 
	  lookup.StartLocsIterator( q_start->Kmer() );
	lookup_table::locs_iterator iloc_end = 
	  lookup.StopLocsIterator( q_start->Kmer() );

	// retrieve all the offsets stored in [iloc_begin,iloc_end) 
	// interval (i.e. all offsets for the current Kmer on the
	// current chunk):

	for ( lookup_table::locs_iterator loc_iter = iloc_begin; 
	      loc_iter != iloc_end; ++loc_iter ) {
	  //	  offsets.push_back( *loc_iter - r );
	  hits.push_back( RawHit( *loc_iter, q_start->QueryPos(), q_start->IsRc() ) );
	}
    } 
}



/// Looks ups queries specified by the iterator interval
/// [q_start, q_end) in the lookup table \c lookup
/// and records all offsets of the \em original query sequences (not queries themselves!). 
/// This method is
/// very similar to QueriesToHits() but is a little bit more economical and faster
/// in cases when the full information stored in hits is not required
/// (such as hit direction, separate offset of the Kmer on the target and offset
/// of the Kmer within the query), but rather a single offset of the original
/// sequence, from which a query was built, on the target suffices.
///
/// NOTE 1: this method retrieves offsets
/// <em>from the currently loaded chunk</em> of the lookup table only!
///
/// NOTE 2: this method \em does \em not check for wrap-arounds: namely,
/// if a query (kmer) passed to this method happens to hit an offset X
/// on the reference, and the offset of that Kmer on the query sequence
/// is Y > X, then the calculated and stored offset of the query sequence on
/// the target will be (unsigned int)(X-Y).
///
/// @param lookup reference genome represented as lookup table
/// @param q_start start of the range of queries
/// @param q_end   end of the range of queries
/// @param[out] hits all hits (occurences on the genome) found for 
/// all queries in the range [q_start, q_end) will be \em appended
/// to this vector (previous content is \em not destroyed).
/// @return total number of hits found.
template < class QueryIterator >
void QueriesToSeqOffsets(const lookup_table & lookup, 
		      QueryIterator q_start,
		      QueryIterator q_end,
		      vec<unsigned int> & offsets) {

 
  for ( ; q_start != q_end ; ++q_start ) {
	  int r = q_start->QueryPos();
     	  unsigned int index = q_start->Kmer();

	  // get the start/stop of locations of this Kmer in the current chunk:
	  lookup_table::locs_iterator iloc_begin = 
	                             lookup.StartLocsIterator(index);
	  lookup_table::locs_iterator iloc_end = 
	                             lookup.StopLocsIterator(index);

	  // retrieve all the offsets stored in [iloc_begin,iloc_end) 
	  // interval (i.e. all offsets for the current Kmer on the
	  // current chunk); shift offsets back by the Kmer position in the 
	  // read (e.g. transform Kmer position on the reference into the
	  // start position of the read alignment suggested by this Kmer):
	  for ( lookup_table::locs_iterator iter = iloc_begin; iter != iloc_end; ++iter ) {
	        // note: the following unsigned int can wrap around if
	        // location offset *iter is less than Kmer offset in query, r
	    	offsets.push_back( *iter - r ); //##############
	  }
  }  // end of loop over all Kmers in the current query sequence

}


#ifdef __GNUC__
// overload of the method for slist
template < class QueryIterator >
void QueriesToSeqOffsets(const lookup_table & lookup, 
		      QueryIterator q_start,
		      QueryIterator q_end,
		      slist<unsigned int> & offsets) {

 
  for ( ; q_start != q_end ; ++q_start ) {
	  int r = q_start->QueryPos();
     	  unsigned int index = q_start->Kmer();

	  // get the start/stop of locations of this Kmer in the current chunk:
	  lookup_table::locs_iterator iloc_begin = 
	                             lookup.StartLocsIterator(index);
	  lookup_table::locs_iterator iloc_end = 
	                             lookup.StopLocsIterator(index);

	  slist<unsigned int>::iterator list_iter = offsets.begin();
	  slist<unsigned int>::iterator list_end = offsets.end();
	  // retrieve all the offsets stored in [iloc_begin,iloc_end) 
	  // interval (i.e. all offsets for the current Kmer on the
	  // current chunk); shift offsets back by the Kmer position in the 
	  // read (e.g. transform Kmer position on the reference into the
	  // start position of the read alignment suggested by this Kmer):
	  for ( lookup_table::locs_iterator iter = iloc_begin; iter != iloc_end; ++iter ) {
	    unsigned int offs = *iter-r;
	    while ( list_iter != list_end ) {
	      if ( *list_iter < offs ) list_iter++;
	      else break;
	    }
	    if ( list_iter != list_end && *list_iter == offs ) continue; // we already have this offset!
	    offsets.insert( list_iter, offs); // we ddid not have this offset, insert it!
	  }
  }  // end of loop over all Kmers in the current query sequence

}
#endif

/// Same as BasesToQueries(const vecbasevector&, VecQueryVec &, unsigned int,
/// AlignDir) (see the description for full details), but takes as its arguments
/// lookup table (Kmer size is taken from that table) and, additionally, 
/// a frequency cutoff \c maxFreq. The output collection of queries will contain
/// only the queries for Kmers with a) non-zero frequencies in the specified 
/// lookup table, and b) frequencies lower than the specified cutoff. If
/// maxFreq=0 (default), no cutoff will be applied (but still only Kmers with non-zero
/// frequencies will be collected). The collection of queries (templatized) can be any
/// 2-d collection of Query objects (VecQueryVec, vec <vector <Query> > , etc) that
/// supports ::resize() and ::operator[] on the first dimension, and ::reserve(unsigned int) 
/// and push_back(Query) on the second dimension.
template < class Container2D >
void BasesToQueries(const vecbasevector &bases, Container2D &queries, 
		    const lookup_table & lookup,
		    unsigned int maxFreq = 0 ,
		    AlignDir direction = FW_OR_RC) {

  int npasses = (( direction == FW ) ? 1 : 2) ;
  unsigned int K = lookup.K();

  queries.resize(bases.size());
  KmerIndexSeq next(K);

  basevector b_tmp;

  for ( int i = 0 ; i < bases.size() ; i++ ) {
    const basevector & b = bases[i];
    if ( b.size() < K ) continue; // no queries would fit, go to next sequence
    const unsigned int last = b.size() - K + 1; // total # of queries on b (one-way)
    queries[i].reserve(npasses*last);
    next.Reset(b);
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      unsigned int index = next(b);
      unsigned int freq = lookup.Freq(index);
      if ( freq == 0 ) continue; // Kmer not found in lookup table, skip it
      if ( maxFreq > 0 && freq >= maxFreq ) continue; // frequency too high, skip it
      queries[i].push_back( Query( j , index ) );
    }
    if ( npasses == 1 ) continue; // only FW, go get next sequence
    // add queries for the rc sequence:
    b_tmp.ReverseComplement(b);
    next.Reset(b_tmp);
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      unsigned int index = next(b_tmp);
      unsigned int freq = lookup.Freq(index);
      if ( freq == 0 ) continue; // Kmer not found in lookup table, skip it
      if ( maxFreq > 0 && freq >= maxFreq ) continue; // frequency too high, skip it
      queries[i].push_back( Query( j , index, true ) ); // add rc query
    }
  }
}



/// Same as BasesToQueries(const vecbasevector &, Container2D &,
/// const lookup_table &, unsigned int, AlignDir), but always assumes that
/// both forward and reverse complement queries are requested (direction=FW_OR_RC)
/// and, for each sequence, stores these queries separately in the corresponding
/// arrays fw_queries and rc_queries. Takes arbitrary 2-d containers of Query objects
/// (such as VecQueryVec, vec<vec<Query> >,
/// etc) that support ::resize() and operator[] on the first dimension, and ::reserve() and 
/// ::push_back() on the second dimension. Non-zero <start_pos> parameter is equivalent to
/// (virtual) truncating the reads: only bases [start_pos, read.size()) will be used to generate
/// both forward and reverse queries and offsets of the kmer queries will be computed \em with 
/// \em respect \em to \em start_pos (i.e. first forward kmer at position start_pos in the read 
/// will have offset 0). The actual sequences in <bases>, however, will not be truncated.
template < class Container2D >
void BasesToQueries(const vecbasevector &bases, Container2D &fw_queries, 
		    Container2D &rc_queries, const lookup_table & lookup,
		    unsigned int maxFreq=0, Bool useReverse = False, unsigned int start_pos=0) {

  unsigned int K = lookup.K();

  // we overestimate the size here, since some queries
  // will be discarded later (those with Kmer frequencies == 0
  // or above the cutoff). 

  fw_queries.resize(bases.size());
  rc_queries.resize(bases.size());
  KmerIndexSeq next(K);
  KmerIndexSeq next_rc(K);

  basevector b_rc;

  for ( size_t i = 0 ; i < bases.size() ; i++ ) {
   
    const basevector & b = bases[i];
    if ( b.size() - start_pos < K ) continue; // no queries would fit, go to next sequence
    const unsigned int last = b.size() - K - start_pos + 1; // total # of queries on b (one-way)
    if ( useReverse ) b_rc.Reverse(b);
    else b_rc.ReverseComplement(b);

    fw_queries[i].reserve(last);
    rc_queries[i].reserve(last);
    next.Reset(b,start_pos);
    next_rc.Reset(b_rc); // here we read from the very start (pos=0) of the reverse
                         // (complement) basevector b_rc, which corresponds to the *end*
                         // of the original basevector, so that no offset is required. The
                         // start_pos argument will be respected because we will retrieve only
                         // <last> queries (Kmers)
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      unsigned int index = next(b);
      unsigned int freq = lookup.Freq(index);
      if ( freq != 0 ) {
	if ( maxFreq == 0 || freq <  maxFreq ) {
	  fw_queries[i].push_back( Query( j , index ) );
	}
      }

      index = next_rc(b_rc);
      freq = lookup.Freq(index);
      if ( freq == 0 ) continue; // Kmer not found in lookup table, skip it
      if ( maxFreq > 0 && freq >= maxFreq ) continue; // frequency too high, skip it
      rc_queries[i].push_back( Query( j , index, true ) ); // add rc query
    }
  }
}


#endif

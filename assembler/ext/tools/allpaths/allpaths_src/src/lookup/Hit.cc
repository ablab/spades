///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "lookup/Hit.h"
#include "lookup/LookupTable.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec< Query, MempoolAllocator<Query> >;
template class OuterVec<QueryVec>;

template class SmallVec< ProcessedHit, MempoolAllocator<ProcessedHit> >;
template class OuterVec<ProcessedHitVec>;


/// Translates DNA sequence <bases> into the collection of K-mers (queries) found 
/// at every position in this sequence; new queries are appended to the <queries> vector;
/// each query stores the position i in the 
/// original sequence (offset) and numerical representation (K-mer index) of K
/// bases starting at that position
void BasesToQueries(const basevector &bases, vec<Query> &queries, unsigned int K, unsigned int offset)
{
  KmerIterator<basevector::const_iterator> KmerSeq(K,bases);

  while(KmerSeq.HasNext()) {
      queries.push_back(Query(offset++,KmerSeq.Next()));
  }
}


/// Transform a collection of basevectors (e.g. reads) into a collection of queries. 
/// Each resulting query stores the Kmer (as numeric index) and the position (offset), 
/// at which this Kmer occurs in the corresponding basevector sequence (\em not
/// the position on the concatenated basevector sequence - see overloaded 
/// implementations!!). In this implementation, \c queries is a vector of vectors -
/// queries[i] stores all the queries derived from bases[i]. For each sequence bases[i],
/// the corresponding collection of queries queries[i] is sorted by 1) direction (fw, rc),
/// and then 2) by Kmer offset position (ascending order) within each direction.
///
/// @param[in] bases vector of sequences (e.g. reads) to be transformed to queries
/// @param[out] queries queries[i] is a vector(-like) object that stores all queries
///   generated for the sequence bases[i]; queries[i] can be empty if no queries could
///   be generated for bases[i] sequence.
/// @param K Kmer length
/// @param direction whether to generate queries only for the specified sequences (FW)
///  or for both specified sequences and their reverse complements (FW_OR_RC); in the 
///  latter case, both fw and rc queries for bases[i] will be stored in queries[i],
///  with the directionality flag set for each query.
void BasesToQueries(const vecbasevector &bases, 
		    VecQueryVec &queries,
		    unsigned int K,
		    AlignDir direction) {

  int npasses = (( direction == FW ) ? 1 : 2) ;

  unsigned int total_queries = 0;
  for ( size_t i = 0 ; i < bases.size() ; i++ ) {
    // how many queries can be generated from the i'th sequence:
    total_queries += npasses*max( (unsigned int)0, bases[i].size() - K + 1 );  
  }
  queries.Reserve(total_queries*sizeof(Query), bases.size());
  queries.resize(bases.size());
  KmerIndexSeq next(K);

  basevector b_tmp;

  for ( size_t i = 0 ; i < bases.size() ; i++ ) {
    const basevector & b = bases[i];
    if ( b.size() < K ) continue; // no queries would fit, go to next sequence
    const unsigned int last = b.size() - K + 1; // total # of queries on b (one-way)
    queries[i].reserve(npasses*last);
    next.Reset(b);
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      queries[i].push_back( Query( j , next(b) ) );
    }
    if ( npasses == 1 ) continue; // only FW, go get next sequence
    // add queries for the rc sequence:
    b_tmp.ReverseComplement(b);
    next.Reset(b_tmp);
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      queries[i].push_back( Query( j , next(b_tmp), true ) ); // add rc query
    }
  }
}


/// Same as BasesToQueries(const vecbasevector &, VecQueryVec &, unsigned int, AlignDir),
/// but always assumes that direction is FW_OR_RC (both forward and reverse) and stores separately 
/// all forward and reverse queries for each sequence in the two passed arrays, fw_queries and 
/// rc_queries, respectively.
void BasesToQueries(const vecbasevector &bases, VecQueryVec &fw_queries, VecQueryVec rc_queries,
		    unsigned int K) {
  unsigned int total_queries = 0;
  for ( size_t i = 0 ; i < bases.size() ; i++ ) {
    // how many queries can be generated from the i'th sequence in one direction:
    total_queries += max( (unsigned int)0, bases[i].size() - K + 1 );  
  }
  fw_queries.Reserve(total_queries*sizeof(Query), bases.size());
  fw_queries.resize(bases.size());
  rc_queries.Reserve(total_queries*sizeof(Query), bases.size());
  rc_queries.resize(bases.size());
  KmerIndexSeq next(K);
  KmerIndexSeq next_rc(K);

  basevector b_rc;

  for ( size_t i = 0 ; i < bases.size() ; i++ ) {
    const basevector & b = bases[i];
    if ( b.size() < K ) continue; // no queries would fit, go to next sequence
    b_rc.ReverseComplement(b);
    const unsigned int last = b.size() - K + 1; // total # of queries on b (one-way)
    fw_queries[i].reserve(last);
    rc_queries[i].reserve(last);
    next.Reset(b);
    next_rc.Reset(b_rc);
    for ( unsigned int j = 0 ; j < last ; ++j ) {
      fw_queries[i].push_back( Query( j , next(b) ) );
      rc_queries[i].push_back( Query( j , next_rc(b_rc), true ) ); // add rc query
    }
  }
}


/// Transform a vecbasevector into queries.  Note that \c queries is a linear
/// vector - query positions are \em not simply offsets of Kmers in each query
/// sequnce. Instead, in order to enable re-assigning queries back to the sequences
/// they originate from, the query positions are
/// derived from the concatenation of all the basevectors [which
/// should therefore be less than 2^31 bases in total length] as (basevector
/// position on the concatenated sequence + offset of the Kmer within the current
/// basevector).
void BasesToQueries(const vecbasevector &bases, vec<Query> &queries, unsigned int K)
{
  size_t v;
  unsigned int q=0;
  for (v=0; v<bases.size(); ++v)
    q += max(0, bases[v].isize() - int(K-1));
  queries.reserve(q);
  KmerIndexSeq next(K);
  unsigned int pos = 0;
  for (v=0; v<bases.size(); ++v) {
    const basevector &b = bases[v];
    if (b.size()>=K) {
      const unsigned int last = b.size()-K+1;
      next.Reset(b); 
      for (unsigned int i=0; i<last; ++i) {
	queries.push_back(Query(pos+i, next(b)));
      }
    }
    pos += b.size();
  }
}


void ClusterHits::operator()(ProcessedHitVec &hits)
{
  sort(hits.begin(), hits.end(), CompareProcessedHitsByRcContigQueryStart());
  bool rc;
  unsigned int contig; 
  longlong pos;
  ProcessedHitVec::iterator it, last, out = hits.begin();
  for (it=hits.begin(); it!=hits.end(); it = last, ++out) {
    // Find range of close-enough hits
    rc = it->IsRc();
    contig = it->TargetContig();
    pos = it->QueryStartOnTarget();
    for (last=it+1; last!=hits.end(); ++last) {
      if (last->IsRc() != rc || last->TargetContig() != contig
	  || last->QueryStartOnTarget() > pos + bw)
	break;
      pos = last->QueryStartOnTarget();
    }
    // Now record this cluster.  We're being slightly tricky here
    // because we're overwriting the input but it's safe because out <= it.
    if (it>out) *out = *it;
    out->SetNHits(totalHits(it, last));
    out->SetBandwidth(pos - it->QueryStartOnTarget());
  }
  // When done, shrink list to valid region and sort by quality 
  hits.resize(distance(hits.begin(), out));
  sort(hits.begin(), hits.end(), CompareProcessedHitsByQuality());
}





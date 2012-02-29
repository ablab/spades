///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOKUP_HIT_H
#define LOOKUP_HIT_H

#include "CoreTools.h"
#include "Basevector.h"
#include "lookup/KmerIndex.h"
#include "lookup/PerfectCount.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

class lookup_table;

/// \file Hit.h

/// \brief A Query encodes a kmer from the query sequence (as numeric index),
/// position (offset) of this kmer in the query sequence (0-based), and the
/// orientation of the query sequence (is this Kmer found at the specified
/// offset in the original query sequence or in its reverse complement). 
///
/// A query does not store any information about the identity of the original
/// sequence it was generated from, so it is responsibility of the client
/// to match sequences to the queries. Query fits
/// in 64 bits and works for any reasonable query (length up to 2^31).
///
struct Query {
  unsigned int k; ///< kmer index
  unsigned int query_pos; ///< offset of the kmer in the sequence (top bit==1 encodes rc queries)
  Query(unsigned int in_query_pos = 0, unsigned int k = 0, bool rc = false) :
    k(k),
    query_pos(rc ? (in_query_pos | TopBit32) : in_query_pos)
  {}
  /// Returns true if the query was generated from rc sequence
  bool IsRc() const { return (query_pos & TopBit32)!=0; }
  /// Returns Kmer index represented by this query
  unsigned int Kmer() const { return k; }
  /// Offset of this query's Kmer in the original sequence or in the
  /// reverse complement of the original sequence (if IsRc() is true)
  int QueryPos() const { return (query_pos & Bits31); }
};

TRIVIALLY_SERIALIZABLE(Query);
typedef SerfVec<Query> QueryVec;
typedef MasterVec<QueryVec> VecQueryVec;

inline ostream & operator<<(ostream &out, const Query &h)
{
  return out << h.QueryPos() << (h.IsRc() ? "/" : "\\") << h.Kmer();
}

#define RETURN_IF_UNEQUAL(expr1, expr2) \
{ if ((expr1)<(expr2)) return true; if ((expr2)<(expr1)) return false; }

struct CompareQueriesByPos {
  bool operator()(const Query &a, const Query &b) {
    return (a.QueryPos() < b.QueryPos());
  }
};

struct CompareQueriesByKmer {
  bool operator()(const Query &a, const Query &b) {
    return (a.Kmer() < b.Kmer());
  }
};

/// Transform a single basevector to queries. Each query in the resulting vector stores
/// the Kmer (as numeric index) and the position, at which this Kmer occurs
/// in the \c bases sequence, incremented by constant amount <offset>. 
/// This method <em>does not</em> generate queries coming
/// from the reverse complement of \c bases.
void BasesToQueries(const basevector &bases, vec<Query> &queries, unsigned int K, unsigned int offset = 0);



/// Same as BasesToQueries(const basevector &bases, vec<Query> &queries, unsigned int K)
/// but takes basevector iterator as its argument instead of basevector itself: this 
/// overload can be used, in particular, with adaptors
template <class ITER>
void BasesToQueries(const ITER & bvec_iter, const ITER & bvec_iter_end, vec<Query> &queries, unsigned int K, unsigned int offset = 0)
{
  KmerIterator<ITER> KmerSeq(K,bvec_iter,bvec_iter_end);

  while ( KmerSeq.HasNext() ) {
    queries.push_back( Query( offset++,KmerSeq.Next() ) );
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
void BasesToQueries(const vecbasevector &bases, VecQueryVec &queries, unsigned int K,
		    AlignDir direction=FW_OR_RC);


/// Same as BasesToQueries(const vecbasevector &, VecQueryVec &, unsigned int, AlignDir),
/// but always assumes that direction is FW_OR_RC (both forward and reverse) and stores separately 
/// all forward and reverse queries for each sequence in the two passed arrays, fw_queries and 
/// rc_queries, respectively.
void BasesToQueries(const vecbasevector &bases, VecQueryVec &fw_queries, VecQueryVec & rc_queries,
		    unsigned int K);


/// Transform a vecbasevector into queries.  Note that \c queries is a linear
/// vector - query positions are \em not simply offsets of Kmers in each query
/// sequnce. Instead, in order to enable re-assigning queries back to the sequences
/// they originate from, the query positions are
/// derived from the concatenation of all the basevectors [which
/// should therefore be less than 2^31 bases in total length] as (basevector
/// position on the concatenated sequence + offset of the Kmer within the current
/// basevector).
void BasesToQueries(const vecbasevector &bases, vec<Query> &queries, unsigned int K);

/// A RawHit encodes a position at which a query kmer matches a target
/// kmer.  This structure is closely related to Query (but does not store
/// additional information that could be used to identify the query a hit
/// corresponds to - it is responsibility of the client to keep this
/// association when needed). The information provided by a RawHit is:
///  - an offset into the reference genome, at which the query kmer is
/// found on the reference; 
///  - offset of the query kmer in the original
/// full query sequence or its reverse complement (see next item); 
///  - orientation of the original query sequence from which the query Kmer 
///    was generated (if IsRc()==true, the query sequence was reverse
///    complemented). 
/// RawHit fits in 64 bits and works for any reasonable query:
/// query length is up to 2^31.  
struct RawHit {
  unsigned int offset; 
  //  int query_pos; ///< negative query_pos values are used to encode rc hits
  unsigned int query_pos; ///< top bit==1 in query_pos values encodes rc hits
  RawHit(unsigned int offset, unsigned int in_query_pos, bool rc ) :
    offset(offset),
    //    query_pos(rc ? -in_query_pos -1 : in_query_pos)
    query_pos ( rc ? (in_query_pos | TopBit32) : in_query_pos )
  { }
  RawHit(unsigned int offset = 0, unsigned int in_query_pos = 0) :
    offset(offset),
    //    query_pos(rc ? -in_query_pos -1 : in_query_pos)
    query_pos ( in_query_pos )
  { }
  //  bool IsRc() const { return query_pos<0; }
  bool IsRc() const { return (query_pos & TopBit32)!=0; }
  unsigned int Offset() const { return offset; }
  //  int QueryPos() const { return (query_pos>=0) ? query_pos : -(1+query_pos); }
  unsigned int QueryPos() const { return query_pos & Bits31; }
  longlong QueryStartOnTarget() const { return longlong(Offset()) - longlong(QueryPos()); }
};

inline ostream & operator<<(ostream &out, const RawHit &h)
{
  return out << h.QueryPos() << (h.IsRc() ? ":" : "=") << h.Offset();
}

struct RawHitIsFw {
  bool operator()(const RawHit &a) {
    return !a.IsRc();
  }
};

struct CompareRawHitsByOffset {
  bool operator()(const RawHit &a, const RawHit &b) {
    return (a.Offset() < b.Offset());
  }
};

struct RawHitsEqualByOffset {
  bool operator()(const RawHit &a, const RawHit &b) {
    return (a.Offset() == b.Offset());
  }
};


struct RawHitOffsetIsBefore {
  RawHit b;
  RawHitOffsetIsBefore(RawHit b) : b(b) { } 
  bool operator()(const RawHit &a) {
    return (a.Offset() < b.Offset());
  }
};

struct CompareRawHitsByQueryStartOffset {
  bool operator()(const RawHit &a, const RawHit &b) {
    return (a.QueryStartOnTarget() < b.QueryStartOnTarget());
  }
};

struct RawHitsEqualByQueryStartOffset {
  bool operator()(const RawHit &a, const RawHit &b) {
    return (a.QueryStartOnTarget() == b.QueryStartOnTarget());
  }
};

struct CompareRawHitsByQueryStart {
  bool operator()(const RawHit &a, const RawHit &b) {
    RETURN_IF_UNEQUAL(a.QueryStartOnTarget(), b.QueryStartOnTarget());
    return a.Offset() < b.Offset();
  }
};


/// A ProcessedHit adds to a RawHit some additional information: how
/// many times did this offset occur?  What contig of target does it
/// correspond to?  What is the range of QueryStartOnTarget values
/// being encoded by this object? Fits in 128 bits.

struct ProcessedHit : public RawHit {
  unsigned short nhits;
  unsigned short bw;
  unsigned int contig;
  ProcessedHit(unsigned int offset = 0, unsigned int contig = 0,
	       int in_query_pos = 0, bool rc = false, unsigned int in_nhits = 0) :
    RawHit(offset, in_query_pos, rc),
    nhits(min(in_nhits, static_cast<unsigned int>(numeric_limits<unsigned short>::max()))),
    bw(0),
    contig(contig)
  { }
  unsigned int TargetContig() const { return contig; }
  unsigned short NHits() const { return nhits; }
  unsigned short Bandwidth() const { return bw; }
  void SetNHits(unsigned int n) { nhits = min(n, static_cast<unsigned int>(numeric_limits<unsigned short>::max())); }
  void SetBandwidth(unsigned int b) { bw = min(b, static_cast<unsigned int>(numeric_limits<unsigned short>::max())); }
};

TRIVIALLY_SERIALIZABLE(ProcessedHit);
typedef SerfVec<ProcessedHit> ProcessedHitVec;
typedef MasterVec<ProcessedHitVec> VecProcessedHitVec;

inline ostream & operator<<(ostream &out, const ProcessedHit &h)
{
  return out << static_cast<RawHit>(h) << "(" << h.TargetContig() << ")["
	     << h.NHits() << "," << h.Bandwidth() << "]";
}

struct CompareProcessedHitsByRcContigQueryStart {
  bool operator()(const ProcessedHit &a, const ProcessedHit &b) {
    RETURN_IF_UNEQUAL(a.IsRc(), b.IsRc());
    RETURN_IF_UNEQUAL(a.TargetContig(), b.TargetContig());
    return (a.QueryStartOnTarget() < b.QueryStartOnTarget());
  }
};

// Note that a standard ascending sort puts better hits at the
// beginning
struct CompareProcessedHitsByQuality {
  bool operator()(const ProcessedHit &a, const ProcessedHit &b) {
    return a.NHits() > b.NHits();
  }
};

template <typename It>
unsigned int totalHits(It first, It last)
{
  unsigned int h = 0;
  for ( ; first!=last; ++first) {
    h += first->NHits();
  }
  return h;
}

// Prepare hits for alignment with gaps by grouping nearby hits and
// forming clusters when the hits are close enough together.  We also
// sort the vector to put the most promising hits first, to improve
// efficiency of rejecting bad alignments.
class ClusterHits {
  int bw;
public:
  ClusterHits(int bw) : bw(bw) { }
  // Use this version to increase bandwidth of ProcessedHits clusters
  void operator()(ProcessedHitVec& hits);
  int bandwidth() const { return bw; }
};



#undef RETURN_IF_UNEQUAL

#endif

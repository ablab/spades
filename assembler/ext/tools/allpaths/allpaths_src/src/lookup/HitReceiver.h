/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef HIT_RECEIVER_H
#define HIT_RECEIVER_H

/// File for classes that fit HitReceiver template in LookupTable::FindHits.
/// \file HitReceiver.h
///
/// 
#include "lookup/LookupTable.h"
#include "lookup/LookAlign.h"
#include "lookup/Hit.h"
#include "pairwise_aligners/SmithWatBandedA.h"


class HitReceiverBase {
protected:  
  /// Lookup table we work with
  lookup_table &look;
  /// Query information
  const vecbasevector &bases;

  /// work vector
  basevector S;
  unsigned int query;
  unsigned int queryLength;

  bool rc;
  /// Whether to change queries
  bool QueryChange(unsigned int queryNumber, bool isRc)
  { return (query!=queryNumber || rc != isRc); }
  /// Change queries
  void SelectQuery(unsigned int queryNumber, bool isRc)
  {
    if (QueryChange(queryNumber, isRc)) {
      if (isRc)
	S.ReverseComplement(bases[queryNumber]);
      else 
	S = bases[queryNumber];
      query = queryNumber;
      rc = isRc;
      queryLength = S.size();
    }
  }

  // only used by the GlobalHitReceiver
  int Bandwidth() const { return -1; }
  bool AmbiguousRead(int query_id) const { return false; }


  HitReceiverBase(lookup_table &look, 
		     const vecbasevector &bases):
     look(look), bases(bases), S(bases[0]), query(0), queryLength(bases[0].size()), rc(false) {}
};

/// Fits HitReceiver template, extends hits into alignments with no indels.
/// \class GlobalUngappedHitReceiver
/// Adapted from UniqueGlobalUngappedAlignCollector.
template <class Collector>
class GlobalUngappedHitReceiver: public HitReceiverBase {
private:  
  /// Where we collect the alignment data. 
  /// This is a policy class that also decides which alignments to 
  /// keep and which to throw out. 
  Collector * coll_;

  /// If gap_ is true, use SmithWatermanBanded to find gaps.
  /// No effort is made to group together different clusters, or even to 
  /// remove identical alignments.
  bool gap_; 

  longlong timesCalled_;///< for debugging

  static const int BANDWIDTH=10;///< for SmithWatBanded
  
public:
   GlobalUngappedHitReceiver(lookup_table &look, 
			     const vecbasevector &bases,
			     Collector * coll,
			     bool gap = false):
     HitReceiverBase(look, bases),  coll_(coll), gap_(gap),  timesCalled_(0)
  { coll_->resize(bases.size()); }

  /// This operator meets HitReceiver interface specified by LookupTable.h
  void operator()(unsigned int queryNumber, unsigned int queryPos, bool isRc,
		  unsigned int targetOffset, unsigned int targetContig,
		  unsigned int numberHits);

  void ChunkDone() { /* nothing to do */ }

  ~GlobalUngappedHitReceiver() { PRINT(timesCalled_); }
};

/// This class accepts lookup-table hits and makes ungapped
/// alignments.  It collects only the unique ones into aligns, while
/// tracking the #errors in the best 2 alignments found in best_errors
/// and second_best_errors.
class UniqueGlobalUngappedHitReceiver : public HitReceiverBase {
private:  
  // Where we collect the alignment data.  
  vec<look_align> &aligns;
  vec<int> &best_errors, &second_best_errors;
  // Parameters
  int best_prox, max_errs;
  // Set up for alignment generation.  All the alignments this class
  // generates will be ungapped.
  look_align la;
  /// Internal use for skipping hits on useless queries
  bool ProvedAmbiguous(unsigned int q)
  { return ( best_errors[q] ==0 && second_best_errors[q] <= best_prox ); }
public:
  UniqueGlobalUngappedHitReceiver(lookup_table &look, const vecbasevector &bases,
				  vec<look_align> &aligns, vec<int> &best_errors,
				  vec<int> &second_best_errors, 
				  int best_prox, int max_errs):
    HitReceiverBase(look, bases),  
    aligns(aligns), best_errors(best_errors),
    second_best_errors(second_best_errors),
    best_prox(best_prox), max_errs(max_errs)
  {
    aligns.resize(bases.size());
    best_errors.resize(bases.size(), infinitely_many);
    second_best_errors.resize(bases.size(), infinitely_many),
    la.nhits = la.indels = 0;
    la.a.SetNblocks(1);
    la.a.SetGap( 0, 0 );
  }
  /// This operator meets HitReceiver interface specified by LookupTable.h
  void operator()(unsigned int queryNumber, unsigned int queryPos, bool isRc,
		  unsigned int targetOffset, unsigned int targetContig,
		  unsigned int numberHits);
  void ChunkDone() { /* nothing to do */ }
};

/// This class accepts lookup-table hits and makes general alignments.
/// It delegates to an alignment collector what to do with them. 
template<class Collector>
class GlobalHitReceiver : public HitReceiverBase {
private:  
  /// Whether to accept only aligns at start, and if so how big a
  /// window at start.  Negative means "do not filter", 0 means
  /// "accept only aligns starting at position 0 on some contig", k>0
  /// means "accept only aligns starting at position <=k on some
  /// contig."
  int atStart;
  /// Tool used to process the hits that come in
  ClusterHits cluster;
  /// Where we collect the alignment data. 
  /// This is a policy class that also decides which alignments to 
  /// keep and which to throw out. 
  Collector * coll;
  // Set up for alignment generation.  
  look_align la;
  // work vectors
  basevector R; // Used to hold a piece of reference for SmithWatBandedA
  ProcessedHitVec hits; // Used to hold the hits before & after clustering
  /// Change queries
  void SelectQuery(unsigned int queryNumber, bool isRc)
  {
    if (QueryChange(queryNumber, isRc)) {
      ProcessHits();
      HitReceiverBase::SelectQuery(queryNumber, isRc);
    }
  }
  /// Process the collected hits on a query we're done with.  Only
  /// used internally--called when QueryChange(), and also at the end
  /// of a chunk.
  void ProcessHits();
  /// Create one alignment from one hit.  Factored out of ProcessHits.
  void ProcessHit(const ProcessedHit h);  
public:
  GlobalHitReceiver(lookup_table &look, const vecbasevector &bases,
		    Collector *coll, int bandwidth, int atStart):
    HitReceiverBase(look, bases),
    atStart(atStart),
    cluster(bandwidth),
    coll(coll)
  {
    coll->resize(bases.size());
    la.nhits = 0;
  }

  /// This operator meets HitReceiver interface specified by LookupTable.h
  void operator()(unsigned int queryNumber, unsigned int queryPos, bool isRc,
		  unsigned int targetOffset, unsigned int targetContig,
		  unsigned int numberHits);
  void ChunkDone() { ProcessHits(); }


  int Bandwidth() const { return cluster.bandwidth(); }
int ErrsAllowed(int query_id) const { return coll->MinErrors(query_id); }
bool IsAmbiguous(int query_id) const { return coll->IsAmbiguous(query_id); }


};


//////////////////////////////////////////////////////////////////////
/// Implementations of template methods
//////////////////////////////////////////////////////////////////////


template<class Collector>
void GlobalUngappedHitReceiver<Collector>::operator()
  (unsigned int q, unsigned int queryPos, bool isRc,
   unsigned int offset, unsigned int contig, unsigned int ) {

  ++timesCalled_;
  if ( offset < queryPos + look.BasesStart() ) return; // Only want aligns within current chunk
  unsigned int q_start = 0, t_start = offset - look.BasesStart() - queryPos;

  unsigned int targetPos = offset - look.ContigStart(contig);
  if (targetPos < queryPos) return; // Only want global aligns of query within contig
  int pos = int(targetPos) - int(queryPos); 

  const basevector &R = look.Bases(); // The bases available in current chunk
  SelectQuery(q, isRc); // get member basevector S set appropriately
  if (t_start + queryLength > R.size()) return; // Only want aligns within current chunk
  if (offset - queryPos + queryLength > look.ContigStop(contig)) return; // Only want global query aligns
  //PRINT4(contig, queryLength, look.BasesStart(), look.ContigStop(contig));
  //PRINT5(queryPos, offset, targetPos, pos, t_start);
  // Compute alignment quality.
  int mismatches = 0;
  for ( unsigned int y = 0; y < queryLength; y++ ) {
    if ( S[q_start + y] != R[ t_start + y ] ) {
      ++mismatches;
    }
  }

  look_align la;
  la.query_id = q;
  la.target_id = contig;
  la.query_length = queryLength;
  la.target_length = look.ContigSize(contig);
  la.rc1 = isRc;

  if (gap_) { //try to find gapped alignments
    unsigned int r0 = offset - look.BasesStart() - queryPos; 
    const basevector &R = look.Bases();
    int errors=0;

    // Create alignment.
    SmithWatBandedA(S, R, -r0, 1 + BANDWIDTH, la.a, errors, 0);
    la.ResetFromAlign(la.a, S, R);
  }
  else {
    /// Generate alignment.  All the alignments this class
    /// generates will be ungapped.
    la.nhits = la.indels = 0;
    la.a.SetNblocks(1);
    la.a.SetGap( 0, 0 );

    // Create alignment.
    la.a.Setpos1(q_start);
    la.a.Setpos2(pos);

    la.a.SetLength( 0, queryLength );
    la.mutations = mismatches;
  }

  // Update. 
  coll_->Insert(la);
  //la.PrintReadableBrief(cout);
  return;
}


template<class Collector>
void GlobalHitReceiver<Collector>::operator()(unsigned int q,
					      unsigned int queryPos, bool isRc,
					      unsigned int offset, unsigned int contig,
					      unsigned int nHits)
{
  // if (0==q) PRINT6(q, queryPos, isRc, offset, contig, nHits);
  if (!coll->AlignsWanted(q))
    return;
  if (offset < queryPos + look.BasesStart()) return; // Align starts before current chunk
  if (offset < queryPos + look.ContigStart(contig)) return; // Align starts before contig
  if (atStart >= 0 && offset > queryPos + look.ContigStart(contig) + atStart + cluster.bandwidth())
    return; // Align starts too late
  SelectQuery(q, isRc); // get member basevector S set appropriately
  if (offset - queryPos + queryLength > look.BasesStop()) return; // Align ends after current chunk
  if (offset - queryPos + queryLength > look.ContigStop(contig)) return; // Align ends after contig

  if ( ! coll->AmAmbiguousCollector() )
    hits.push_back(ProcessedHit(offset, contig, queryPos, isRc, nHits));
  else
  {
    ProcessHit(ProcessedHit(offset, contig, queryPos, isRc, nHits));
    coll->Insert(la);
  }
}

inline int Extend(const basevector & S, const basevector & R, int offset) {
  int mismatches=0;
  for ( unsigned int i = 0; i < S.size(); ++i ) {
    if ( S[i] != R[i-offset] ) {
      ++mismatches;
    }
  }
  return mismatches;
}

/// Create the alignment implied by a ProcessedHit.  Works mostly
/// with class members, for speed.
template<class Collector>
void GlobalHitReceiver<Collector>::ProcessHit(const ProcessedHit h)
{
  // Create alignment.
  la.query_id = query;
  la.target_id = h.TargetContig();
  la.query_length = queryLength;
  la.target_length = look.ContigSize(la.target_id);
  la.rc1 = h.IsRc();

  const unsigned int contigStart = look.ContigStart(la.target_id);
  const unsigned int contigStop = look.ContigStop(la.target_id);
 
  //start of alignment.
  unsigned int r0 = h.Offset();
  if (r0 > static_cast<unsigned int>(h.QueryPos()))
    r0 -= h.QueryPos();
  const unsigned int bw = max(1+cluster.bandwidth(),
			      static_cast<int>(h.Bandwidth()));
  //We carefully cut out a small piece of the reference centered on
  //the hit so as to not put too much of a load on SetToSubOf.  The
  //allowable reference is bounded by three constraints: the desired
  //interval r0 + [-queryLength - bw, 3*queryLength + 2*bw), the available
  //bases [look.BasesStart(), look.BasesStop()), and the current
  //contig [contigStart, contigStop). 
  const unsigned int subStop = min(contigStop, look.BasesStop());
  unsigned int subStart = max(contigStart, look.BasesStart());
  unsigned int subLength = queryLength * 3 + 2*bw;
  if (r0 > queryLength + bw) {
    subStart = max(subStart, r0 - queryLength - bw);
  }
  if (subStart + subLength > subStop)
    subLength = subStop - subStart; 

  R.SetToSubOf(look.Bases(), subStart - look.BasesStart(), subLength);
    
  int errors;
  SmithWatBandedA(S, R, -(r0-subStart), bw, la.a, errors, 0);
  la.ResetFromAlign(la.a, S, R);
  la.a.Setpos2(la.a.pos2() + subStart - contigStart);
}





template<class Collector>
void GlobalHitReceiver<Collector>::ProcessHits()
{
  //Note that all coordinates, for now, are from the beginning of the
  //concatenated reference. That is, they are in the same coordinate
  //system as h.Offset().
  if (hits.empty())
    return;
  if (!coll->AlignsWanted(query)) {
    hits.clear();
    return;
  }

  cluster(hits);
  for (unsigned int i=0; i<hits.size(); ++i) { // Walk through clustered hits
    ProcessHit(hits[i]); // --> result in la

    if (la.pos1() > 0 || la.Pos1() < int(queryLength))
      continue; // not a global alignment
    if (atStart>=0 && la.pos2() > atStart)
      continue; // starts too late into contig

    coll->Insert(la);       
    // if (query==0) la.PrintReadableBrief(cout);
  }
  hits.clear();
}

#endif

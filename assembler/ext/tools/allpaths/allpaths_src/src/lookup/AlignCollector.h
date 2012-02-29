/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIGN_COLLECTOR_H
#define ALIGN_COLLECTOR_H

/// File for AlignCollector-interface classes used by HitReceivers.
/// \file AlignCollector.h
///
/// 
#include "lookup/LookupTable.h"
#include "lookup/LookAlign.h"
#include "lookup/HitReceiver.h"
#include "FastIfstream.h"

/**
   Type concept: AlignCollector

   A collector of alignments: for each query it stores the
   (possibly multiple) alignments of that query to the target.
   For example, when storing alignments of reads to the reference
   genome, each query is a read and each target is a <genome part>.

   An AlignCollector basically represents a vec< vec< look_align > >,
   storing a vector of alignments for each query; but also has
   methods for defining what alignments it will accept.
   ImperfectLookup() finds imperfect alignments of queries to target,
   and stores them in the specified AlignCollector.

   Some classes that model AlignCollector are: UniqueByErrDiffAlignCollector,
   MaxErrDiffAlignCollector.

   A class that is a model of AlignCollector must define the following methods:

   Method: clear

   Erases all alignments from this collector
   
   > void clear();

   Method: resize

   Sets the size of the collector to i and initializes "best errors" to
   "infinitely many" (only for newly added elements!).
   
   > void resize(unsigned int i);

   Method: size
   
   Total size of the collector container (NOT the number of 
   "good" alignments received or all alignments received!). 
   See AlignedCount() and UniquelyAlignedCount().
   
   >unsigned int size() const;

   Method: AlignsWanted
   
   Whether to skip further alignments for this query.
   Reimplemented by unique align collectors.
   
   >bool AlignsWanted(unsigned int query_id) const;

   Method: ErrorThreshold
   
   Returns the maximum number of errors
   an alignment should have in order to be accepted
   or at least considered by the collector. It is
   guaranteed that an alignment with more mismatches than
   the number returned by this method will be rejected.
   
   > unsigned int ErrorThreshold(unsigned int query_id);

   Method: Aligned

   Returns true if an alignment with maxErrs or fewer 
   errors was seen by the collector for the specified query.
   The alignment does not have to be unique and there is even no guarantee
   that the alignment is actually stored in the collector
   (depends on the policy implemented for the specific collector type).
   
   > bool Aligned(unsigned int query_id) const;

   Method: Consolidate
   
   Do nothing for base implementation.
   Some align collectors may need to do some cleanup before
   the results are usable (e.g. sort and/or remove duplicates)
   
   > void Consolidate();

   Method: CollectorName

   Returns the name of this AlignCollector. 

   > const char *GetCollectorName();

*/

/// This class encapsulates the fundamentals of saving alignments:
/// track number of queries and best # errors per query, and
/// parameters regarding the number of allowed errors.
///
/// Abstract base class that helps implement classes that model
/// the type concept AlignCollector.
///
/// \class AlignCollectorBase
///
// TODO: potentially dangerous truncation of index in this class
// and all derived classes due to query_id being passed as int or uint
class AlignCollectorBase {

public:
  
  AlignCollectorBase(int maxErrDiff, int maxErrs, unsigned int maxReads = 0): 
    maxErrDiff_(maxErrDiff), maxErrs_(maxErrs), betterBest_(0)
  { if (maxReads > 0) { this->resize(maxReads); } } 

  virtual ~AlignCollectorBase() {}

  /// Sets the size of the collector to i and initializes "best errors" to
  /// "infinitely many" (only for newly added elements!).
  virtual void resize(unsigned int i) { best_.resize(i, infinitely_many); }

  /// Erases all alignments from this collector
  virtual void clear() { best_.clear(); }
  
  /// Total size of the collector container (NOT the number of 
  /// "good" alignments received or all alignments received!).
  /// See AlignedCount() and UniquelyAlignedCount().
  /// In other words: this returns the number of queries for
  /// which this AlignCollector is storing alignments; there might be
  /// no alignments stored for some queries and multiple alignments
  /// stored for others, depending on the concrete type of the align collector
  /// used and the parameters it was initialized with.
  virtual unsigned int size() const { return best_.size(); }
  
  /// Whether we want to see further alignments for this query.
  /// Returns 'false' to say we do not want to see any more alignments
  /// for this query (i.e. the concrete collector type is designed to collect
  /// unique alignments, and from the history of alignments it has seen so far
  /// it can be guaranteed that no unique alignment can exist for the given query)
  /// Reimplemented by unique align collectors.
  virtual bool AlignsWanted(unsigned int query_id) const { return true; }

  /// Returns the maximum number of errors
  /// an alignment should have in order to be accepted
  /// or at least considered by the collector. It is
  /// guaranteed that an alignment with more mismatches than
  /// the number returned by this method will be rejected.
  virtual unsigned int ErrorThreshold(unsigned int query_id) {
     return infinitely_many;
  }

  /// Returns true if an alignment with maxErrs or fewer 
  /// errors was seen by the collector for the specified query.
  /// The alignment does not have to be unique and there is even no guarantee
  /// that the alignment is actually stored in the collector
  /// (depends on the policy implemented for the specific collector type). 
  virtual bool Aligned(unsigned int query_id) const {
    return best_[query_id] <= maxErrs_;
  }
  
  /// Queries if there is a unique alignment for query <query_id> in the collector;
  /// it is left for the collector to decide what is unique or whether to reimplement this method at all.
  virtual bool UniquelyAligned(unsigned int query_id) const {
      return false;
  }

  /// returns an align for query <id>, if such align exists; behavior if alignment does
  /// not exist is left to concrete implementations; if multiple alignments for the query
  /// are stored by concrete collector implementation, it is up to the collector to decide
  /// which one to return, but it is normally expected to be the best one - see concrete implementations
  /// (some implementations may not be able to return the best alignment correctly until they are consolidated).
  virtual const look_align & Align(unsigned int id) const = 0;

  /// returns mutable align reference that can be changed from external code; use wisely, it's generally a bad idea
  /// to change aligns right inside the collector.
  virtual look_align & MutableAlign(unsigned int id) const = 0;

  /// Print parseable and readable (if <readable>=<true>) alignments for query_id.
  virtual void Print(ostream & out, int query_id, bool readable = true ) const = 0;

  /// Print parseable and readable brief (if <readable> = <true>)  alignments for all query_ids.
  /// Result depends on concrete collector type implementing the method: some collectors may store
  /// and print only unique alignments, others would print multiple alignments per read etc.
  virtual void Print(ostream & out, bool readable = true) = 0;

  /// Do nothing for base implementation.
  /// Some align collectors may need to do some cleanup before
  /// the results are usable (e.g. sort and/or remove duplicates)
  virtual void Consolidate() {}

  // just for completeness
  virtual bool AmAmbiguousCollector() const { return false; }
  virtual bool IsAmbiguous(int query_id) const { return false; }
  virtual void PrintAmb(ostream & out, int query_id, int offset) {};
  virtual void PrintAmb(ostream & out, int offset) {};

  /// pass alignment to the collector; depending on the collector's policy it may decide
  /// to accept and store the alignment or to ignore it.
  virtual void Insert(const look_align & la) = 0;

  // Delete all alignments associated with <query_id>.
  virtual void Erase(int query_id) = 0;

  /// Adjust the query_ids by adding offset to query_id of each alignment stored.
  virtual void AdjustQueryIds(int offset) = 0;

  /// Adjust the target_ids by adding offset to target_id of each alignment stored.
  virtual void AdjustTargetIds(int offset) = 0;

  /// Number of errors in best alignment for query  <query_id> or -1 if no alignment was ever received.
  virtual int MinErrors(int query_id) const { 
    return (best_[query_id] < (unsigned int)infinitely_many ? (int)best_[query_id] : -1); 
  }

  /// Print list of MinErrors values for all query_ids, in minAlignErrors.txt format.
  virtual void PrintMinErrors(ostream & out) {
    for (int i=0; i != best_.isize(); ++i ) out << MinErrors(i) << "\n";
  }

  static const char *GetCollectorName() { return "AlignCollectorBase"; }

protected:
  /// Protected field: maxErrDiff_
  /// We will reject alignments that have more than maxErrDiff_ more
  /// errors than the best alignment.  (Note that we may also
  /// reject the best alignment if it has more than maxErrs_ errors.)
  unsigned int maxErrDiff_;
  
  /// Protected field: maxErrs_
  /// We will reject alignments with more than this many errors --
  /// as if the alignment was not there.
  unsigned int maxErrs_;
  
  /// Protected field: best_
  /// For each query, the number of errors in the best
  /// alignment of that query to the target.
  vec<unsigned int> best_;
  
  longlong betterBest_;

};  // class AlignCollectorBase


/// This class encapsulates saving only aligns that are good enough (multiple
/// alignment can be "good enough" and thus stored by this collector!).
/// \class MaxErrDiffAlignCollector
///
/// An alignment is good enough for this collector if
///    * it is the best possible alignment and it has best_errs <= maxErrs errors, OR
///    * it is not the best alignment, but it has errs <= best_errs + maxErrDiff
/// 
/// This is one example of a policy class that also decides 
/// which alignments to keep and which to throw out. 
///
/// A different implementation could use an ordered list
/// internally for each query instead of a vec, which would probably save
/// memory and eliminate the need to Consolidate().
class MaxErrDiffAlignCollector : public AlignCollectorBase {

public:
  typedef vec<look_align>::const_iterator ConstIter;
  typedef vec<look_align>::iterator Iterator;
  
  MaxErrDiffAlignCollector(int maxErrDiff, int maxErrs, unsigned int maxReads = 0): 
    AlignCollectorBase(maxErrDiff, maxErrs), consolidated_(true) // empty collector is consolidated!
    { if ( maxReads > 0 ) resize(maxReads); } 

  /// Sets the size of the collector to i and initializes "best errors" to
  /// "infinitely many" (only for newly added elements!).
  void resize(unsigned int i) { 
    aligns_.resize(i); 
    AlignCollectorBase::resize(i); 
  }
  
  /// Erases all elements of the collector container
  void clear() { aligns_.clear(); AlignCollectorBase::clear(); }

  /// Add in a look_align (or don't if we don't want it).
  void Insert(const look_align & la);

  /// Returns the maximum number of errors
  /// an alignment should have in order to be accepted
  /// or at least considered by the collector. It is
  /// guaranteed that an alignment with more mismatches than
  /// the number returned by this method will be rejected.
  unsigned int ErrorThreshold(unsigned int query_id) {
      return best_[query_id]+maxErrDiff_;
  }

  /// Clean up: remove any look_aligns that are not good enough.  Sort
  /// the remaining ones.  Idempotent: if called again, does nothing.
  void Consolidate();

  /// Beginning of the look_align container for this query_id
  ConstIter Begin(int query_id) const { 
    ForceAssert(consolidated_);
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  ConstIter End(int query_id) const { 
    ForceAssert(consolidated_);
    return aligns_[query_id].end(); 
  }

  /// Beginning of the look_align container for this query_id
  Iterator Begin(int query_id) { 
    ForceAssert(consolidated_);
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  Iterator End(int query_id)  { 
    ForceAssert(consolidated_);
    return aligns_[query_id].end(); 
  }

  /// Number of alignments for that query_id
  int Size(int query_id) const { 
    ForceAssert(consolidated_);
    return aligns_[query_id].size(); 
  }

  /// Returns immutable collection of aligns stored for the
  /// specified query id.
  const vec<look_align> & Aligns(int query_id) const {
    ForceAssert(consolidated_);
    return aligns_[query_id];
  }

  /// If the collector is consolidated and there is still at least one
  /// align for the query_id, then returns the best align; otherwise breaks.
  const look_align & Align(unsigned int query_id ) const {
    ForceAssert(consolidated_);
    ForceAssertGt(aligns_[query_id].size(),0u);
    return aligns_[query_id][0];
  }

  /// best alignment for the query that can be modified (use wisely!); 
  /// if the query did not have any alignments stored, this method will crash on ForceAssert 
  look_align & MutableAlign(unsigned int query_id) const { 
    ForceAssert(consolidated_);
    ForceAssertGt(aligns_[query_id].size(),0u);
    return const_cast<look_align &>(aligns_[query_id][0]);
  }

  /// Adjust the query_ids by adding offset to query_id of each alignment stored.
  void AdjustQueryIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustQueryId( i, offset);
    }
  }

  /// Adjust the target_ids by adding offset to target_id of each alignment stored.
  void AdjustTargetIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustTargetId( i, offset);
    }
  }


  /// Returns true if a unique alignment with n_err <= maxErr
  /// errors was seen by this collector for the specified query.
  /// An alignment is considered to be unique 
  /// if no other alignment for the same query with n_err+maxErrDiff
  /// or fewer errors was seen.  
  bool UniquelyAligned(unsigned int query_id) const {
    // this collector always keeps all multiple alignments
    // with err <= best_err+maxErrDif, so that the alignment
    // is unique (so far) only if the number of stored alignments
    // for the query == 1; the first check below allows using this
    // method with non-consolidated collectors - those are guaranteed
    // to keep only alignments below best_err + maxErrDiff cutoff, but
    // are not guaranteed that best <= maxErr is satisfied.
    return best_[query_id] <= maxErrs_ && aligns_[query_id].size() == 1;
  }

  /// Print all parseable and readable brief (if <readable>=<true>) alignments 
  /// up to best_err+maxErrDiff for 
  /// read <query_id> if the read aligns (best alignment has best_err <= maxErr errors);
  /// otherwise print nothing at all.
  void Print(ostream & out, int query_id, bool readable = true ) const {
    if ( ! Aligned(query_id ) ) return;
    for (ConstIter i=Begin(query_id); i != End(query_id); ++i) {
	i->PrintParseable(out);
	if ( readable ) i->PrintReadableBrief(out);
    }
  }

  /// Print all parseable and readable brief (if <readable> = <true> )
  /// alignments up to best_err+maxErrDiff for all query_ids that 
  /// align (have best align errors best_err <= maxErrs).
  /// Calls Consolidate() first (just in case).
  void Print(ostream & out, bool readable = true) {
    Consolidate();
    for (int i=0; i != aligns_.isize(); ++i ) Print(out, i, readable);
  }

  // Delete all alignments associated with query_id.
  void Erase(int query_id);

  static const char *GetCollectorName() { return "MaxErrDiffAlignCollector"; }
  

private:
  bool consolidated_;
  vec<vec<look_align> > aligns_;

  /// Clean up: remove look_aligns for that query that are not good enough.
  void EraseBad(int query_id);

  /// Adjust query id by adding offset to alignment at specified position.
  void AdjustQueryId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->query_id += offset;
    }
  }


  /// Adjust target id by adding offset to alignments at specified position.
  void AdjustTargetId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->target_id += offset;
    }
  }

  
  
};  // class MaxErrDiffAlignCollector

/// This class encapsulates saving only aligns for queries that align uniquely.
/// \class UniqueByErrDiffAlignCollector
///
/// By uniquely we mean that the next best alignment has at least
/// maxErrDiff+1 more errors than in the best possible alignment, and
/// also that the best alignment has at most maxErrs errors. 
/// 
class UniqueByErrDiffAlignCollector : public AlignCollectorBase {

public:
  typedef vec<look_align>::const_iterator Iterator;
  
  UniqueByErrDiffAlignCollector(int maxErrDiff, int maxErrs, unsigned int maxReads = 0): 
    AlignCollectorBase(maxErrDiff, maxErrs)
    { if ( maxReads > 0 ) resize(maxReads); } 
   

  /// Sets the size of the collector to i and initializes all "errors" to
  /// "infinitely many" (only for newly added elements!).
  void resize(unsigned int i) { 
      aligns_.resize(i); 
      second_best_.resize(i,infinitely_many); 
      AlignCollectorBase::resize(i); 
    }

  /// Erases all the elements of this collector container.
  void clear() {
    aligns_.clear();
    second_best_.clear();
    AlignCollectorBase::clear();
  }

  // Erases an entry of id i
  void Erase(int i) {
          best_[i] = infinitely_many;
          second_best_[i] = infinitely_many;
  }
  
  /// Whether to accept more alignments for this query: returns false
  /// if the query is already proved ambiguous. The query is deemed ambiguous
  /// if a) the current alignment is ambiguous (the difference between
  /// numbers of mismatches in the best and next best alignments found so far
  /// is maxErr or less) \em and b) the best alignment found so far has
  /// maxErr mismatches or less [it is easy to see that under these two 
  /// conditions, no new unique alignment can be found]
  bool AlignsWanted(unsigned int query_id) const
  { 
    return ! ( best_[query_id] <= maxErrDiff_ &&
	       second_best_[query_id] - best_[query_id] <= maxErrDiff_ );
  }

  /// Add in a look_align (or don't if we don't want it).
  void Insert(const look_align & la);

  /// Returns the maximum number of errors
  /// an alignment should have in order to be accepted
  /// or at least considered by the collector. It is
  /// guaranteed that an alignment with more mismatches than
  /// the number returned by this method will be totally rejected and will
  /// have no side effects (e.g. such as updating the second-best error rate).
  /// Note that the alignment with ErrorThreshold() or 
  /// fewer errors might still not be accepted.
  unsigned int ErrorThreshold(unsigned int query_id) {
    return NextErrors(query_id) - 1 ;
  }

  /// Whether this query_id has a unique alignment with <= maxErrs errors
  bool UniquelyAligned(unsigned int query_id) const { 
    return best_[query_id] <= maxErrs_
      && second_best_[query_id] > maxErrDiff_ + best_[query_id]; 
  }

  /// Number of alignments for this query_id is 0 or 1
  int Size(int query_id) const { return UniquelyAligned(query_id); }
 
  /// best alignment for the query; if the query did not have a unique
  /// alignment, the returned value may contain junk (but the call is safe).
  const look_align & Align(unsigned int query_id) const { 
    return aligns_[query_id];
  }

  /// best alignment for the query that can be modified (use wisely!); 
  /// if the query did not have a unique
  /// alignment, the returned value may contain junk (but the call is safe).
  look_align & MutableAlign(unsigned int query_id) const { 
    return const_cast<look_align &>(aligns_[query_id]);
  }

  /// Adjust the query_ids by adding offset to each one.
  void AdjustQueryIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      aligns_[i].query_id += offset;
    }
  }


  /// Adjust the target_ids by adding offset to each one.
  void AdjustTargetIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      aligns_[i].target_id += offset;
    }
  }

  /// Returns the number of errors in the second best alignment
  /// found so far for the query sequence query_id.
  int NextErrors(int query_id) const {
    return second_best_[query_id];
  }

  /// Print parseable and readable brief (if <readable> = <true> ) unique alignment for query_id
  /// if such alignment exists; otherwise prints nothing at all.
  void Print(ostream & out, int query_id, bool readable = true) const {
    if (UniquelyAligned(query_id)) {
      aligns_[query_id].PrintParseable(out);
      if ( readable ) aligns_[query_id].PrintReadableBrief(out);
    }
  }

  /// Print parseable and readable brief (if <readable> = <true>) unique alignments for all query_ids
  /// that align uniquely.
  void Print(ostream & out, bool readable = true) {
	for (int i=0; i != aligns_.isize(); ++i ) Print(out, i, readable);
  }

  static const char *GetCollectorName() { return "UniqueByErrDiffAlignCollector"; }

private:
  vec<look_align> aligns_;
  vec<unsigned int> second_best_;
  
};  // class UniqueByErrDiffAlignCollector





/// This policy encapsulates saving only the best and next best aligns as long as best has 
/// maxErr or less errors. The "second best" is defined as the best addition(s) to the *single* best 
/// alignment, i.e. if there are multiple "best" alignments with the same number of errors,
/// all but one will be considered "second best" (the "best" assignment will be made based
/// on the order collector receives the alignments), and no alignments with greater numbers
/// of errors will be stored. If there is indeed a single best alignment, then all alignments
/// with next-best number of errors will be stored as "second best".
///
/// This class differs from UniqueByErrDiffAlignCollector in that 1) it actually stores
/// the second best alignment(s) rather than only mismatch count for them; 2) it does not
/// enforce any errDiff filtering and does not implement early refusal to accept more alignments,
/// so it is guaranteed to store the true best and next best
/// alignments as long as best err <= maxErr.
/// [ Note that UniqueByErrDiffAlignCollector does refuse to accept more
/// alignments when it has enough information to decide that the read can not
/// be aligned uniquely. Hence UniqueByErrDiffAlignCollector's "best" align
/// stored for non-uniquely aligning reads can be incorrect.]
///
/// This collector always keeps aligns for a given read sorted: best alignment first,
/// followed by next best alignment(s). The Align() method returns the best alignment
/// (compatible with UniqueByErrDiff...).
///
/// While the collection logic does not use ErrDiff, the 
/// \em definitions of aligned and uniquely aligned reads are the same as those for
/// UniqueByErrDiffAlignCollector. Namely, the accessor logical
/// methods Aligned(read_id) and UniquelyAligned(read_id) are still available and
/// return true under exactly the same conditions (see UniqueByErrDiffAlignCollector).
/// \class BestNextBestAlignCollector
///

class BestNextBestAlignCollector : public AlignCollectorBase {

public:
  typedef vec<look_align>::const_iterator ConstIter;
  
  BestNextBestAlignCollector(int maxErrDiff, int maxErrs, unsigned int maxReads = 0): 
    AlignCollectorBase(maxErrDiff, maxErrs)
    { if ( maxReads > 0 ) resize(maxReads); } 
 

  /// Sets the size of the collector to i and initializes all "errors" to
  /// "infinitely many" (only for newly added elements!).
  void resize(unsigned int i)
  { aligns_.resize(i); second_best_.resize(i,infinitely_many); 
    AlignCollectorBase::resize(i); }

  /// Erases all the elements of this collector container.
  void clear() {
    aligns_.clear();
    second_best_.clear();
    AlignCollectorBase::clear();
  }
  
  /// Add in a look_align (or don't if we don't want it).
  void Insert(const look_align & la);


  /// Returns the maximum number of errors
  /// an alignment should have in order to be accepted
  /// or at least considered by the collector. It is
  /// guaranteed that an alignment with more mismatches than
  /// the number returned by this method will be rejected.
  /// For this specific policy implementation, an alignment
  /// with ErrorThreshold() or fewer errors is indeed guaranteed
  /// to be accepted.
  unsigned int ErrorThreshold(unsigned int query_id) {
    return NextErrors(query_id) ;
  }

  /// Clean up: removes duplicates (only second best alignments can be 
  /// duplicated in this collector without cleanup)
  void Consolidate();


  /// Whether this query_id has a unique alignment
  bool UniquelyAligned(unsigned int query_id) const { 
    return best_[query_id] <= maxErrs_
      && second_best_[query_id] > maxErrDiff_ + best_[query_id]; 
  }

  /// Number of alignments for this query_id is 0 or 1
  int Size(int query_id) const { return aligns_[query_id].isize(); }
 
  /// Best alignment for this query_id. WARNING: if no alignments at
  /// all are stored for this query id the method will crash on ForceAssert!
  /// Check for existence of at least one alignment prior to invoking this method.
  const look_align & Align(unsigned int query_id) const { 
    ForceAssertGt(aligns_[query_id].size(),0u);
    return aligns_[query_id][0];
  }

  /// best alignment for the query that can be modified (use wisely!); 
  /// if the query did not have any alignments stored, this method will crash on ForceAssert 
  look_align & MutableAlign(unsigned int query_id) const { 
    ForceAssertGt(aligns_[query_id].size(),0u);
    return const_cast<look_align &>(aligns_[query_id][0]);
  }

  /// Beginning of the look_align container for this query_id
  ConstIter Begin(int query_id) const { 
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  ConstIter End(int query_id) const { 
    return aligns_[query_id].end(); 
  }

  /// Adjust the query_ids by adding offset to each one.
  void AdjustQueryIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustQueryId( i, offset);
    }
  }

  /// Adjust the target_ids by adding offset to each one.
  void AdjustTargetIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustTargetId( i, offset);
    }
  }

  /// Returns the number of errors in the second best alignment
  /// found so far for the query sequence query_id.
  int NextErrors(int query_id) const {
    return second_best_[query_id];
  }

  /// Print parseable and readable brief (if <readable>=<true>) best and next best 
  /// alignments for query_id if this query_id aligns (has best_ett <= maxErr).
  void Print(ostream & out, int query_id, bool readable = true) const {
    if ( ! Aligned( query_id ) ) return;
    for (ConstIter i=Begin(query_id); i != End(query_id); ++i) {
	i->PrintParseable(out);
	if ( readable ) i->PrintReadableBrief(out);
    }
  }

  /// Print parseable and readable brief (if <readable>=<true>) best and next best 
  /// alignments for all query_ids that align (have best_err <= maxErr).
  void Print(ostream & out, bool readable = true) {
      Consolidate();
      for (int i=0; i != aligns_.isize(); ++i ) Print(out, i, readable);
  }

  // Delete all alignments associated with query_id.
  void Erase(int q) {
    best_[q] = infinitely_many;
    second_best_[q] = infinitely_many;
    aligns_[q].clear();
  }

  static const char *GetCollectorName() { return "BestNextBestAlignCollector"; }

private:
  vec< vec<look_align> > aligns_;
  vec<unsigned int> second_best_;
  
  typedef vec<look_align>::iterator Iterator;

  /// Beginning of the look_align container for this query_id
  Iterator Begin(int query_id) { 
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  Iterator End(int query_id)  { 
    return aligns_[query_id].end(); 
  }

  /// Adjust query id by adding offset to alignment at the given position
  void AdjustQueryId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->query_id += offset;
    }
  }

  /// Adjust target id by adding offset to alignment at the given position.
  void AdjustTargetId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->target_id += offset;
    }
  }
};  // class BestNextBestAlignCollector




/// Always keeps the best align(s) (by error count) as long as they have <= maxErr; the 
/// distance maxErrDiff is ignored by this policy except for the AlignedUniquely() method 
/// (this collector still holds the error count for next-best align observed, so it can tell
/// whether the best align it holds is unique or not). The difference from UniqueByErrDiffAlignCollector
/// is that truly best alignment with <= maxErr errors (or multiple best alignments if there are many) observed so far 
/// is guaranteed to be held regardless of its uniqueness status.
class BestAlignCollector : public AlignCollectorBase {

public:
  typedef vec<look_align>::const_iterator ConstIter;
  
  BestAlignCollector(int maxErrDiff, int maxErrs, unsigned int maxReads = 0): 
    AlignCollectorBase(maxErrDiff, maxErrs)
    { if ( maxReads > 0 ) resize(maxReads); } 
 

  /// Sets the size of the collector to i and initializes all "errors" to
  /// "infinitely many" (only for newly added elements!).
  void resize(unsigned int i)
  { aligns_.resize(i); second_best_.resize(i,infinitely_many); 
    AlignCollectorBase::resize(i); }

  /// Erases all the elements of this collector container.
  void clear() {
    aligns_.clear();
    second_best_.clear();
    AlignCollectorBase::clear();
  }
  
  /// Add in a look_align (or don't if we don't want it).
  void Insert(const look_align & la);


  /// Returns the maximum number of errors
  /// an alignment should have in order to be accepted
  /// or at least considered by the collector. It is
  /// guaranteed that an alignment with more mismatches than
  /// the number returned by this method will be rejected.
  /// For this specific policy implementation, an alignment
  /// with ErrorThreshold() or fewer errors is indeed guaranteed
  /// to be accepted.
  unsigned int ErrorThreshold(unsigned int query_id) {
    return NextErrors(query_id) ;
  }

  /// Clean up: removes duplicates 
  void Consolidate();


  /// Whether this query_id has a unique alignment
  bool UniquelyAligned(unsigned int query_id) const { 
    return best_[query_id] <= maxErrs_
      && second_best_[query_id] > maxErrDiff_ + best_[query_id]; 
  }

  /// Number of alignments for this query_id is 0 or 1
  /// (DJ: This isn't true - a number larger than one can be returned.)
  int Size(int query_id) const { return aligns_[query_id].isize(); }
 
  /// Best alignment for this query_id. WARNING: if no alignments at
  /// all are stored for this query id the method will crash on ForceAssert!
  /// Check for existence of at least one alignment prior to invoking this method.
  /// (DJ: I'm not sure what this returns if there is more than one best align.)
  const look_align & Align(unsigned int query_id) const { 
    ForceAssertGt(aligns_[query_id].size(),0u);
    return aligns_[query_id][0];
  }

  /// best alignment for the query that can be modified (use wisely!); 
  /// if the query did not have any alignments stored, this method will crash on ForceAssert 
  look_align & MutableAlign(unsigned int query_id) const { 
    ForceAssertGt(aligns_[query_id].size(),0u);
    return const_cast<look_align &>(aligns_[query_id][0]);
  }

  /// Beginning of the look_align container for this query_id
  ConstIter Begin(int query_id) const { 
    return aligns_[query_id].begin(); 
  }
  ConstIter Begin2(int query_id) const { 
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  ConstIter End(int query_id) const { 
    return aligns_[query_id].end(); 
  }
  ConstIter End2(int query_id) const { 
    return aligns_[query_id].end(); 
  }

  /// Adjust the query_ids by adding offset to each one.
  void AdjustQueryIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustQueryId( i, offset);
    }
  }

  /// Adjust the target_ids by adding offset to each one.
  void AdjustTargetIds(int offset) {
    for (unsigned int i=0; i != size(); ++i) {
      AdjustTargetId( i, offset);
    }
  }

  /// Returns the number of errors in the second best alignment
  /// found so far for the query sequence query_id.
  int NextErrors(int query_id) const {
    return second_best_[query_id];
  }

  /// Print all parseable and readable brief (if <readable>==<true>) best alignments for 
  /// the query_id if the read aligns (best_err <=maxErr, not necesserily uniquely), or does nothing
  /// if read does not align/aligns with > maxErr errors.
  void Print(ostream & out, int query_id, bool readable = true) const {
    if ( ! Aligned(query_id) ) return;
    for (ConstIter i=Begin(query_id); i != End(query_id); ++i) {
	i->PrintParseable(out);
	if ( readable ) i->PrintReadableBrief(out);
    }
  }

  /// Print all parseable and readable brief (if <readable> = <true>) best alignments for all query_ids
  /// that align (i.e. have at least one alignment with err <= maxErr)
  void Print(ostream & out, bool readable = true) {
      Consolidate();
      for (int i=0; i != aligns_.isize(); ++i ) Print(out, i, readable);
  }

  // Delete all alignments associated with query_id.
  void Erase(int q) {
    best_[q] = infinitely_many;
    second_best_[q] = infinitely_many;
    aligns_[q].clear();
  }

  static const char *GetCollectorName() { return "BestAlignCollector"; }

private:
  vec< vec<look_align> > aligns_;
  vec<unsigned int> second_best_;
  
  typedef vec<look_align>::iterator Iterator;

  /// Beginning of the look_align container for this query_id
  Iterator Begin(int query_id) { 
    return aligns_[query_id].begin(); 
  }

  /// End of the look_align container for this query_id.
  Iterator End(int query_id)  { 
    return aligns_[query_id].end(); 
  }

  /// Adjust query id by adding offset to alignment at the given position
  void AdjustQueryId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->query_id += offset;
    }
  }

  /// Adjust target id by adding offset to alignment at the given position.
  void AdjustTargetId(int id, int offset) {
    for (Iterator i=Begin(id); i != End(id); ++i) {
      i->target_id += offset;
    }
  }
};  // class BestAlignCollector



// This class compares alignments of a query sequence to a known alignment for that query
// and decides if the query is "ambiguous" or not.
// \class AmbiguousAlignCollector
//
// Takes as input a limit L for each read [equal to the # errors previously found in an alignment 
// to the target], output a boolean value for each read: is it ambiguous or unambiguous?  
// 
// A read is deemed ambiguous if this collector receives either (a) one alignment of the same
// read with <L errors; or (b) two alignments of the read with L errors.  
// 
// Alignments of the read with >L errors are irrelevant to the output.

class  AmbiguousAlignCollector : public AlignCollectorBase
{
public:

AmbiguousAlignCollector(int maxErrDiff, int maxErrs, vec<int> & thrsPerRead)
  :  AlignCollectorBase(maxErrDiff, maxErrs, 0)
{
  counts_per_read_.clear();
  counts_per_read_.resize(thrsPerRead.size(),0);
  num_better_aligns_=2;
  best_.resize(thrsPerRead.size());
  for ( unsigned int i = 0 ; i < thrsPerRead.size() ; i++ ) {
    best_[i] = thrsPerRead[i]; // reuse best_ to save thresholds as initial best errors
  }
} 

/// Ambiguous collector can not be resized; only i==current size is allowed
void resize(unsigned int i)
{
  ForceAssertEq( i , size() );
}

/// not implemented for this collector; will print error message and crash.
  const look_align & Align(unsigned int i) const { cout << "AmbiguousAlignCollector::Align() is not implemented" << endl; ForceAssert(0); 
    look_align* la = new look_align;
    return *la;    }
/// not implemented for this collector; will print error message and crash.
look_align & MutableAlign(unsigned int i) const { cout << "AmbiguousAlignCollector::MutableAlign() is not implemented" << endl; ForceAssert(0); 
    look_align* la = new look_align;
    return *la;    }
/// not implemented for this collector; won't do anything
void Print(ostream & s, int id, bool readable ) const {}; 
/// not implemented for this collector; won't do anything
void Print(ostream & s, bool readable ) {}; 
/// not implemented for this collector; will print error message and crash.
void Erase(int id) { cout << "AmbiguousAlignCollector::Erase() is not implemented" << endl; ForceAssert(0); }
/// not implemented for this collector; won't do anything
void AdjustQueryIds(int offset) {};
/// not implemented for this collector; won't do anything
void AdjustTargetIds(int offset) {};


/// Erases all elements of this collector (size will be 0!)
void clear() {
   AlignCollectorBase::clear();
   best_.clear();
}

bool AmAmbiguousCollector() const { return true; }


/// Returns the maximum number of errors
/// an alignment should have in order to be accepted
/// or at least considered by the collector. It is
/// guaranteed that an alignment with more mismatches than
/// the number returned by this method will be rejected.
/// For this specific policy implementation, an alignment
/// with ErrorThreshold() or fewer errors is not guaranteed
/// to be accepted (it will not if the read is already
/// deemed ambiguous), AlignsWanted() takes precedence in this decision.
unsigned int ErrorThreshold(unsigned int query_id) {
   return best_[query_id] ;
}

// Use align to test read for ambiguity, add align if needed
void Insert(const look_align & la);

bool AlignsWanted(unsigned int query_id) const { 
  return (counts_per_read_[query_id] < num_better_aligns_); 
}
bool IsAmbiguous(int query_id) const { 
  return  (counts_per_read_[query_id] ==  num_better_aligns_); 
}

void SetAsAmbiguous(int query_id) { counts_per_read_[query_id] = num_better_aligns_; }

int MinErrors(int query_id) const { return best_[query_id]; }

void AddAmbiguousReads(vec<int> &amb_reads, int first_read, int last_read) {
  for (int i=first_read; i < last_read; ++i )
  {
    if ( counts_per_read_[i-first_read] == num_better_aligns_ )
      amb_reads.push_back(i);
  }
}

void PrintMinErrors(ostream & out) {
  for (int i=0; i != best_.isize(); ++i ) {
    out << best_[i] << "\n";  
  }
}

// Print ambiguous for query_id
void PrintAmb(ostream & out, int query_id, int offset) {
  if ( IsAmbiguous(query_id) ) {
    int q(query_id+offset);
    out << q << endl;
  }
}

// Print ambiguous for all query_ids
void PrintAmb(ostream & out, int offset) {
  for (int i=0; i != counts_per_read_.isize(); ++i ) PrintAmb(out, i, offset);
}

  static const char *GetCollectorName() { return "AmbiguousAlignCollector"; }
  
private:
// ordered by query_id
// vec<int> threshold_per_read_;
 vec<int> counts_per_read_;
int num_better_aligns_;
  
};  // class AmbiguousAlignCollector


/// Load look_aligns into a look_align collector: specific collector class
/// instance defines what aligns are acceptable and have to be kept. Only alignments 
/// for the reads with Query ID (as recorded in the qltout file) in the range 
/// [<low_id>, <high_id>) will be loaded with defaults <low_id>=0, <high_id>=max(QueryID).
/// NOTE: if <low_id> is non-zero, the Query ID of the loaded alignments will be adjusted
/// down by -<low_id> (i.e. if <low_id>=1000, the alignment for query with ID=1000 will
/// be recorded in the collector as having query id=0, use ::AdjustQuerytIds() if needed).
/// The last parameter tells whether alignments in the qltout file are ordered by query_id
/// (default false). This option has effect only if \c high_id is set and allows to abort reading
/// large qltout files as soon as alignment for query \c high_id is reached.
/// Collector must be resized prior to passing to this method, so that its size is 
/// (at least) equal to the number
/// of query IDs in the qltout file or in the [<low_id>,<high_id>) range (if specified)! 
/// Upon return, the align collector
/// is consolidated and ready to use.
template <class AlignCollector>
void LoadLookAligns(const String & file_name, AlignCollector &aligns, int low_id = 0, int high_id = -1, bool ordered=false ) {
  String line;
  fast_ifstream in(file_name);
  bool good = false;
  look_align la;
  unsigned int count = 0;
  const String tag("QUERY");
  while(1) {
    good = getline_if_match( in, line , tag);
    if ( in.fail( ) ) break;
    if ( good ) {
      la.ReadParseable(line);
      if ( low_id > 0 && la.QueryId() < low_id ) continue;
      if ( high_id >=0 && la.QueryId() >= high_id ) {
	if ( ordered ) break;
	continue;
      }
      la.query_id -= low_id;
      aligns.Insert(la);
      count++;
    }
  }
  aligns.Consolidate();
}

/// Same as LoadLookAligns() (see docs) but reads from a binary qltb file
template <class AlignCollector>
void LoadLookAlignsBinary(const String & file_name, 
			  AlignCollector &aligns, 
			  int low_id = 0, 
			  int high_id = -1, 
			  bool ordered=false ) {

  vec<look_align> a;

  // if ordered, will use these for searching:
  look_align low;
  low.query_id=low_id;
  look_align high;
  high.query_id = high_id-1;
  
  longlong size = BinarySize3<look_align>(file_name, true);

  longlong chunk = 1000000; // million at a time

  bool done = false;

  for ( longlong from = 0 ; from < size && ! done ; from += chunk ) {

    longlong to = min ( from+chunk, size );
    BinaryReadRange3(file_name,from,to,a,true);
    
    vec<look_align>::iterator it;
    vec<look_align>::iterator last;
    if ( low_id > 0 && ordered ) {
      it = lower_bound(a.begin(), a.end(), low, order_lookalign_Query() );
      if ( it == a.end() ) continue; // all query ids in this chunk are < low_id
    } else it = a.begin();
    if ( high_id >= 0 && ordered ) {
      last = upper_bound(a.begin(), a.end(), high, order_lookalign_Query() );
      if ( last != a.end() ) done = true; // highest query id in this chunk is >= high_id
    } else last = a.end();

    for ( ; it != last ; ++it ) {
      if ( low_id > 0  && it->QueryId() < low_id ) continue;
      if ( high_id >=0 && it->QueryId() >= high_id ) continue;
      it->query_id -= low_id;
      it -> PrintReadableBrief(cout);
      aligns.Insert(*it);
    }
  }
  aligns.Consolidate();
}


/// Same as LoadLookAligns() except that it is assumed that alignments are saved in terms of the "local"
/// reference (i.e. amplicons or hybrid selection targets). <l2gmap> provides the mapping from these amplicons
/// onto the whole genome: <l2gmap[i].first> and <l2gmap[i].second> are zero-based contig and offset, respectively, 
/// of the start of the contig (amplicon) found at position <i> in the local reference. Each loaded alignment
/// is translated into the whole-genome coordinates using this map before attempting to add it to the collector.
/// NOTE: this may have non-trivial effect if the amplicons of the local reference overlap and non-unique
/// alignments against the local reference are stored: a read from the region covered by overlapping
/// amplicons will align to both amplicons (i.e. it is "non-unique", quite artificially, with respect  
/// to the local reference), but when mapped to the global reference these alignments will be equivalent (map 
/// to the same genomic position). If the alignments to the local reference are loaded "as is" and Unique...
/// align collector is used, then such alignments will be lost (align collector will deem them non-unique
/// and discard them). However, if the aligns are loaded with mapping, the collector will be aware of the fact
/// that it is actually the same align.
/// The aligns loaded by this method into the collector will have their TargetId() and StartOnTarget()
/// method returning position in terms of the whole genome. NOTE: Target length will be set to 0 in the 
/// remapped aligns!!!
template <class AlignCollector>
void LoadMappedLookAligns(const String & file_name, vec<pair<unsigned int, unsigned int> > l2g_map, 
			  AlignCollector &aligns, int low_id = 0, int high_id = -1, bool ordered=false ) {
  String line;
  fast_ifstream in(file_name);
  bool good = false;
  look_align la;
  unsigned int count = 0;
  const String tag("QUERY");
  while(1) {
    good = getline_if_match( in, line , tag);
    if ( in.fail( ) ) break;
    if ( good ) {
      la.ReadParseable(line);
      if ( low_id > 0 && la.QueryId() < low_id ) continue;
      if ( high_id >=0 && la.QueryId() >= high_id ) {
	if ( ordered ) break;
	continue;
      }
      la.query_id -= low_id;
      // remap:
      la.SetStartOnTarget( l2g_map[la.target_id].second + la.StartOnTarget() );
      la.target_id = l2g_map[la.target_id].first;
      la.target_length = 0; // THIS IS UNSAFE! CALLER MUST UPDATE TARGET_LENGTH IF NEEDED
      // insert remapped alignment; collector should recognize the attempt to insert the 
      // same align again and should not consider this as inserting a "non-unique" alignment
      aligns.Insert(la);
      count++;
    }
  }
  aligns.Consolidate();
}


/// Same as LoadMappedLookAligns() (see docs) but reads from a binary qltb file
template <class AlignCollector>
void LoadMappedLookAlignsBinary(const String & file_name, 
				vec<pair<unsigned int, unsigned int> > l2g_map, 
				AlignCollector &aligns, 
				int low_id = 0, 
				int high_id = -1, 
				bool ordered=false ) {

  vec<look_align> a;

  // if ordered, will use these for searching:
  look_align low;
  low.query_id=low_id;
  look_align high;
  high.query_id = high_id-1;
  
  longlong size = BinarySize3<look_align>(file_name, true);

  longlong chunk = 1000000; // million at a time

  bool done = false;

  for ( longlong from = 0 ; from < size && ! done ; from += chunk ) {

    longlong to = min ( from+chunk, size );
    BinaryReadRange3(file_name,from,to,a,true);
    
    vec<look_align>::iterator it;
    vec<look_align>::iterator last;
    if ( low_id > 0 && ordered ) {
      it = lower_bound(a.begin(), a.end(), low, order_lookalign_Query() );
      if ( it == a.end() ) continue; // all query ids in this chunk are < low_id
    } else it = a.begin();
    if ( high_id >= 0 && ordered ) {
      last = upper_bound(a.begin(), a.end(), high, order_lookalign_Query() );
      if ( last != a.end() ) done = true; // highest query id in this chunk is >= high_id
    } else last = a.end();

    for ( ; it != last ; ++it ) {
      if ( low_id > 0  && it->QueryId() < low_id ) continue;
      if ( high_id >=0 && it->QueryId() >= high_id ) continue;
      it->query_id -= low_id;
      // remap:
      it->SetStartOnTarget( l2g_map[it->TargetId()].second + it->StartOnTarget() );
      it->target_id = l2g_map[it->TargetId()].first;
      it->target_length = 0; // THIS IS UNSAFE! CALLER MUST UPDATE TARGET_LENGTH IF NEEDED
      // insert remapped alignment; collector should recognize the attempt to insert the 
      // same align again and should not consider this as inserting a "non-unique" alignment
      aligns.Insert(*it);
    }
  }
  aligns.Consolidate();
}



/// Utility method that scans the align collector and returns the number
/// of aligned queries. Normally a query is considered to be aligned if at least
/// one alignment with MaxErr or fewer errors was observed, uniqueness
/// is not taken into account - unless specific collector
/// implementation overrides its Aligned(int) method to give a different
/// definition of being "aligned".
template <class AlignCollector>
unsigned int AlignedCount(const AlignCollector & aligns) {
  unsigned int aligned = 0;
  for ( unsigned int id = 0 ; id < aligns.size() ; id++ ) {
    if ( aligns.Aligned(id) ) aligned++;
  }
  return aligned;
}

/// Utility method, scans the align collector and returns the number
/// of uniquely aligned queries. Normally, a query is considered to be 
/// uniquely aligned if its best alignment has MaxErr or fewer errors, 
/// and if the next best alignment has more than MaxErr+ErrDiff errors -
/// unless specific collector implementation overrides its UniquelyAligned
/// method to give a different definition of being aligned.
template <class AlignCollector>
unsigned int UniquelyAlignedCount(const AlignCollector & aligns) {
  unsigned int uniquely_aligned = 0;
  for ( unsigned int id = 0 ; id < aligns.size() ; id++ ) {
    if ( aligns.UniquelyAligned(id) ) uniquely_aligned++;
  }
  return uniquely_aligned;
}

/// Combines alignments stored in two collectors into one (<trg>);
/// collector's own logic is used to decide whether to accept
/// an alignment stored in <src> into <trg> etc. Collectors must
/// have same size and aligns at same position (read id) from both
/// collectors will be combined and the winner will be stored at the
/// same position in trg

template <class AlignCollector>
void Combine(AlignCollector & trg, const AlignCollector & src ) {
  if ( trg.size() != src.size() ) {
    cout << "Sizes of align collectors do not match" << endl;
    exit(1);
  }
  for ( unsigned int i = 0 ; i < trg.size() ; i++ ) {
    if ( ! trg.AlignsWanted(i) ) continue; // we know all we wanted to know and do not need more alignments
    if ( ! src.Aligned(i) ) continue; // read was not aligned in src at all, skip it
    trg.Insert(src.Align(i)); // ok, try to insert best alignment from src
  }
}

#endif

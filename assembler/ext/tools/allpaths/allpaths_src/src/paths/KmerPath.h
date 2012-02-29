/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef KMER_PATH_H
#define KMER_PATH_H

#include "CoreTools.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "math/Functions.h"
#include "String.h"
#include "paths/KmerPathInterval.h"
#include "CommonSemanticTypes.h"
#include "HashSimple.h"
#include "graph/Digraph.h"

class KmerPathLoc;  // forward declaration

/**
   Class: KmerPath

   A path of kmers where each kmer is shifted by one base relative to the 
   previous kmer.  For example (using 4-mers):

   (begin example)
   AGTC
    GTCT
     TCTA
       ....
   (end example)

   The path is represented as a sequence of <kmer numbers>, rather than the 
   actual sequence of the path;  use <KmerBaseBroker> to get the sequence.   
   The sequence of kmer numbers is represented using ranges:

      >[23-30][55-67][4-4][13-15] ...

   Each range is represented as one <KmerPathInterval>.  A KmerPathInterval 
   may also represent a gap, though right now (4-03-07) this functionality is 
   not used in ALLPATHS.

   To iterate over the kmers in a KmerPath without worrying where <segments> 
   start and end, use <Begin()>/<End()>;  see <KmerPathLoc>.

   There may be several ways to represent the same kmer path in this 
   representation; but see <Canonicalize()>.

   See also: <vecKmerPath>, <KmerPathInterval>
*/
class KmerPath : public SerfVec<KmerPathInterval> {
    typedef MempoolAllocator<KmerPathInterval> Alloc;
    typedef SerfVec<KmerPathInterval> Base;

public:

  typedef KmerPathInterval value_type;

  KmerPath( )
  : SerfVec<KmerPathInterval>()
  { }

  KmerPath( String s );

  KmerPath(const vec<kmer_id_t>& path );

  KmerPath( Alloc alloc )
  : Base(alloc) {}

  int NSegments( ) const { return size( ); }
  Bool InRange( int i ) const { return (0<=i && i<NSegments()); }

  void Clear( ) { clear( ); }

  bool IsEmpty( ) const { return empty( ); }

  KmerPathLoc BasePosToSegmentPos( int base_pos ) const;

  void SetNSegments( int n ) { resize(n); }
  void Reserve( int n ) { reserve(n); }

  const KmerPathInterval& Segment( int i ) const { return (*this)[i]; }

  const KmerPathInterval& FirstSegment( ) const { return (*this)[0]; }
  const KmerPathInterval& LastSegment( ) const { return (*this)[ (int) size( ) - 1 ]; }

  Bool isGap( int i ) const { return (InRange(i) && (*this)[i].isGap( )); }
  Bool isSeq( int i ) const { return (InRange(i) && (*this)[i].isSeq( )); }

  void Reverse( );

  // Method: GetKmer
  // Get the nth kmer from a gap-free KmerPath:
  kmer_id_t GetKmer( int n ) const;

  // Method: Canonicalize
  // Merge consecutive gaps and consecutive mergable sequence segments.
  void Canonicalize( );

  kmer_id_t Start( int i ) const { return (*this)[i].Start( ); }
  kmer_id_t Stop( int i ) const { return (*this)[i].Stop( ); }
  int Length( int i ) const { return (*this)[i].Length( ); }
  int TotalLength( ) const {
    int len = 0;
    for (int ii=0; ii<this->NSegments( ); ii++) len += this->Length( ii );
    return len;
  }

  kmer_id_t Start() const { return FirstSegment().Start(); }
  kmer_id_t Stop() const { return LastSegment().Stop(); }

  // same functions with different names, for thinking about gaps:
  kmer_id_t Minimum( int i ) const { return (*this)[i].Start( ); }
  kmer_id_t Maximum( int i ) const { return (*this)[i].Stop( ); }
  int Stretch( int i ) const { return (*this)[i].Stretch( ); }

  // AddSegment will concatenate abutting KmerPathIntervals.
  void AddSegment( KmerPathInterval rpi );
  void AddSegmentNoConcatenate( KmerPathInterval rpi );
  void AddSegment( kmer_id_t start, kmer_id_t stop, Bool is_gap = False )
  { AddSegment( KmerPathInterval(start, stop, is_gap) ); }
  void AddGap( kmer_id_t start, kmer_id_t stop)
  { AddSegment( KmerPathInterval(start, stop, True) ); }

  void Append( const KmerPath& other, int begin=0, int end=-1 );
  void append( const KmerPath& other, int begin=0, int end=-1 )
  {    return Append( other, begin, end );    }
  void AppendNoFirstKmer( const KmerPath& other, int begin=0, int end=-1 );
  // Don't call AppendNoFirstKmer on a KmerPath starting with a gap!

  void SetStart( int i, kmer_id_t x )
    { (*this)[i].Set( x, Stop(i), isGap(i) ); }
  void SetStop( int i, kmer_id_t x )
    { (*this)[i].Set( Start(i), x, isGap(i) ); }
  // same functions with different names, for thinking about gaps:
  void SetMinimum( int i, kmer_id_t x )
    { (*this)[i].Set( x, Maximum(i), isGap(i) ); }
  void SetMaximum( int i, kmer_id_t x )
    { (*this)[i].Set( Minimum(i), x, isGap(i) ); }
  void SetGap( int i, kmer_id_t x1, kmer_id_t x2 )
    { (*this)[i].Set( x1, x2, True ); }


  // Group: Things mentioning KmerPathLocs.
  //
  // Actual code can't appear until after
  // the complete declaration of KmerPathLoc, so Subpath()s are in
  // KmerPath.cc, and the others appear below (for inlining).

  // NOTE: End() points to the last k-mer, not one-past-the-end.
  // Two KmerPathLoc's represent a closed interval, not a half-open one.
  // (Hence the need for the ...NoFirstKmer variant.  Oh well.)
  KmerPathLoc Begin() const;
  KmerPathLoc End() const;

  // Method: CopySubpath
  // appends subpath to answer; loc1 and loc2 are locs in *this
  void CopySubpath( KmerPathLoc loc1, KmerPathLoc loc2, KmerPath& ans ) const;
  void CopySubpathNoFirstKmer( KmerPathLoc loc1, KmerPathLoc loc2,
			       KmerPath& ans ) const;
  void CopySubpathNoLastKmer( KmerPathLoc loc1, KmerPathLoc loc2,
			      KmerPath& ans ) const;

  void CopyHead( const KmerPathLoc& loc, KmerPath& ans ) const;
  void CopyTail( const KmerPathLoc& loc, KmerPath& ans ) const;
  void CopyHeadNoLastKmer( const KmerPathLoc& loc, KmerPath& ans ) const;
  void CopyTailNoFirstKmer( const KmerPathLoc& loc, KmerPath& ans ) const;
  // Copy the interval [loc1,loc2] adjusting the individual gap sizes
  // so that it fits in the space of a gap with bounds given_{min,max}
  void CopySubpathAdjustGaps( KmerPathLoc loc1, KmerPathLoc loc2,
			      int given_min, int given_max,
			      KmerPath& ans ) const;

  // Method: KmerCount
  // Return the number of kmers in this kmer path, ignoring any gaps.
  // renamed from Length(), to avoid confusion: this ignores gaps.
  int KmerCount( ) const {
    int sum = 0;
    for ( int i = 0; i < NSegments( ); i++ )
      if (isSeq(i)) sum += Length(i);
    return sum;
  }

  int MinLength(int i) const { return ( isSeq(i) ? Length(i) : Minimum(i) ); }
  int MaxLength(int i) const { return ( isSeq(i) ? Length(i) : Maximum(i) ); }

  int AveLength(int i) const { return ( MinLength(i) + MaxLength(i) ) / 2; }

  int MinLength() const { return MinLength(0, NSegments()-1); }
  int MaxLength() const { return MaxLength(0, NSegments()-1); }

  // Length of the path if all the gaps are as small as they can be
  int MinLength(int seg1, int seg2) const {
    int sum = 0;
    for ( int i = seg1; i <= seg2; ++i )
      sum += MinLength(i);
    return sum;
  }

  // Length of the path if all the gaps are as large as they can be
  int MaxLength(int seg1, int seg2) const {
    int sum = 0;
    for ( int i = seg1; i <= seg2; ++i )
      sum += MaxLength(i);
    return sum;
  }

  float MidLength( ) const
  {    float len = 0.0;
       for ( int u = 0; u < NSegments( ); u++ )
       {    const KmerPathInterval& x = Segment(u);
            if ( x.isSeq( ) ) len += x.Length( );
            else len += float( x.Maximum( ) + x.Minimum( ) ) / 2.0;    }
       return len;    }

  Bool GapFree( ) const
  {    for ( int i = 0; i < NSegments( ); i++ )
            if ( isGap(i) ) return False;
       return True;    }

  // Predicate: Proper
  // A KmerPath is "proper" if it is nonempty, has no gaps on either end, and it
  // does not have adjacent segments which are mergeable.
  Bool Proper( ) const;

  // TODO: potentially dangerous truncation of index by next three methods
  void AppendToDatabase( vec<tagged_rpint>& segs, int i ) const;
  void AppendToDatabase( vec<big_tagged_rpint>& segs, int i ) const;
  void AppendToDatabase( vec<new_tagged_rpint>& segs, int i ) const;

  // Method: IsSubpathAnchoredLeft
  // Determine if path 1 formally matches subset of path 2, anchored on left:
  // >   ------1-------
  // >   ----------2----------
  // Read code to see exact handling of subtleties.
  friend Bool IsSubpathAnchoredLeft( const KmerPath& p1, const KmerPath& p2 );

  friend void BinaryWrite( int fd, const KmerPath& p );
  friend void BinaryRead( int fd, KmerPath& p );

  // Pretty-printed output:
  friend ostream& operator<<(ostream& out, const KmerPath& p);
  friend void PrintFolded( ostream& out, const KmerPath& p );

  #define KMER_PATH_HASH_CONST 123456789
  ulonglong GetHash(ulonglong seed=0) const
  {
    for ( int i = 0; i < this->NSegments(); ++i )
      seed = seed * KMER_PATH_HASH_CONST + this->Segment(i).GetHash();
    return seed;
  }

  friend Bool operator==( const KmerPath& p1, const KmerPath& p2 ) {
    if( p1.NSegments( ) != p2.NSegments( ) ) return False;
    for( int i=0; i<p1.NSegments( ); i++ )
      if( ! (p1.Segment(i)==p2.Segment(i)) ) return False;
    return True;
  }

  friend Bool operator!=( const KmerPath& p1, const KmerPath& p2 ) {
    return !(p1==p2);
  }

  // Sorting according to this will leave == KmerPaths adjacent:
  friend Bool operator<( const KmerPath& p1, const KmerPath& p2 ) {
    if( p1.NSegments( ) != p2.NSegments( ) )
      return( p1.NSegments( ) < p2.NSegments( ) );
    for( int i=0; i<p1.NSegments( ); i++ )
      if( ! (p1.Segment(i)==p2.Segment(i)) )
	return( p1.Segment(i)<p2.Segment(i) );
    // If we get here, these are == (not <).
    return False;
  }

  friend Bool operator<=( const KmerPath& p1, const KmerPath& p2 ) {
    if( p1.NSegments( ) != p2.NSegments( ) )
      return( p1.NSegments( ) < p2.NSegments( ) );
    for( int i=0; i<p1.NSegments( ); i++ )
      if( ! (p1.Segment(i)==p2.Segment(i)) )
	return( p1.Segment(i)<p2.Segment(i) );
    // If we get here, these are == (not <).
    return True;
  }
};

SELF_SERIALIZABLE(KmerPath);

inline nbases_t Kmers2Bases( nkmers_t nkmers, nbases_t K ) {
  return nkmers + K - 1;
}

inline nkmers_t Bases2Kmers( nbases_t nbases, nbases_t K ) {
  return nbases - K + 1;
}

/*
   Type: vecKmerPath

   A vector of <kmer paths>.  Can represent, for example, all unipaths computed from a set of reads
   for a given kmer size.  Note that this representation is in <kmer space>;
   <vecbasevector> represents a vector of basevectors in <sequence space>.
*/
typedef MasterVec<KmerPath> vecKmerPath;

// Function: VecOfKmerPath
// Convert a <vecKmerPath> to an ordinary vec<KmerPath>.
inline
vec<KmerPath> VecOfKmerPath( const vecKmerPath& v )
{    vec<KmerPath> p;
     p.resize( v.size( ) );
     for ( size_t i = 0; i < v.size( ); i++ )
       p[i] = v[i];
     return p;    }

// Function: FeudalOfKmerPath
// Convert a vec<KmerPath> to a vecKmerPath.
inline
vecKmerPath FeudalOfKmerPath( const vec<KmerPath>& v )
{
  vecKmerPath paths;
  size_t nseqs = v.size();
  longlong rawsize = 0;
  for (size_t ii=0; ii<nseqs; ii++)
    rawsize += v[ii].NSegments( );
  paths.Reserve( rawsize + nseqs, nseqs );
  for (size_t ii=0; ii<v.size(); ii++)
    paths.push_back( v[ii] );
  return paths;
}

// Returns a vector contain lengths of the corresponding base sequences

inline
vec<int> KmerPathSeqLength( const vecKmerPath& v, const int K )
{
  vec<int> lengths(v.size());
  for (size_t i = 0; i < v.size(); i++)
    lengths[i] = v[i].TotalLength() + (K - 1);
  return lengths;
}

/**
    Class: KmerPathLoc

    Holds a location in a KmerPath, iterator-like.  Philosophically,
    some methods treat it as pointing at one segment of the path,
    others as pointing at a particular kmer in that segment.
    Note that unlike with STL iterators, a range of locations
    [from,to] specified by two KmerPathLoc's includes *both* its
    endpoints.

    KmerPathLoc supports moving forward and backward in the path
    by given amounts ( see IncrementHaltAtGap() and operator+() );
    unlike STL's random access iterator, these are not quite
    constant-time as they depend on the number of KmerPathIntervals
    skipped by the increment/decrement operation.  But with a good
    kmer numbering, the number of these intervals would hopefully
    be relatively small.

    Actually like a const_iterator, and probably I should have
    two versions, a non-const one also.  For now, if you need to
    recover from a KmerPathLoc a non-const KmerPath, you'll have
    to explicitly cast the return value of GetPath() or GetPathPtr().

    Some methods point into gaps in the following way:
    If isGap(), then m_loc holds an integer which indicates how many gap kmers
    you have already passed:
      if m_loc >= 0, then this is measured from the left of the interval;
      if m_loc < 0 then -m_loc-1 is number of gap-mers passed on the right.
    So m_loc=0 means the left edge of a gap, m_loc=-1 means right edge.

    This class is delicate: if you use an in-sequence method while you're
    in a gap, it may silently give you garbage answers.
*/
class KmerPathLoc {
 public:

  KmerPathLoc(const KmerPath& path, int index, int loc=0)
    : mp_path(&path), m_index(index), m_loc(loc) { }
  KmerPathLoc(const KmerPath* pathp, int index, int loc=0)
    : mp_path(pathp), m_index(index), m_loc(loc) { }
  KmerPathLoc(const KmerPathLoc& other)
    : mp_path(other.mp_path), m_index(other.m_index), m_loc(other.m_loc) { }
  KmerPathLoc( ) : mp_path(0), m_index(0), m_loc(0) { }

  void Set(const KmerPath& path, int index, int loc)
  {    mp_path = &path;
       m_index = index;
       m_loc = loc;    }

  const KmerPath& GetPath() const { return *mp_path; }
  const KmerPath* GetPathPtr() const { return mp_path; }
  int GetIndex() const { return m_index; }
  int GetLoc() const { return m_loc; }
  pair<int,int> Where() const { return make_pair(GetIndex(),GetLoc()); }
  KmerPathInterval GetSegment(int off=0) const
  { return mp_path->Segment(m_index+off); }
  // Are you on the first/last segment in the KmerPath?
  Bool atFirst() const { return (m_index == 0); }
  Bool atLast() const { return (m_index == mp_path->NSegments() - 1); }
  // Are you at the left/right edge of your interval?
  Bool atIntervalLeft() const { return (m_loc==0); }
  Bool atIntervalRight() const { return (m_loc==Length()-1); }
  // Are you at the first/last base of the entire KmerPath?
  Bool atBegin() const { return (atFirst() && atIntervalLeft()); }
  Bool atEnd() const { return (atLast() && atIntervalRight()); }
  // Is there a gap just before/after this kmer?
  Bool GapToLeft() const { return (atIntervalLeft() && isGap(-1)); }
  Bool GapToRight() const { return (atIntervalRight() && isGap(+1)); }

  Bool isSeq(int off=0) const { return mp_path->isSeq(m_index+off); }
  Bool isGap(int off=0) const { return mp_path->isGap(m_index+off); }
  // If you're in sequence, use these:
  kmer_id_t Start(int off=0) const { return mp_path->Start(m_index+off); }
  kmer_id_t Stop(int off=0) const { return mp_path->Stop(m_index+off); }
  longlong Length(int off=0) const { return mp_path->Length(m_index+off); }
  kmer_id_t GetKmer() const { return(mp_path->Start(m_index)+m_loc); }
  // If you're in a gap, use these:
  kmer_id_t Minimum(int off=0) const { return mp_path->Minimum(m_index+off); }
  kmer_id_t Maximum(int off=0) const { return mp_path->Maximum(m_index+off); }
  // Use these in either situation:
  longlong MinLength(int off=0) const { return mp_path->MinLength(m_index+off); }
  longlong MaxLength(int off=0) const { return mp_path->MaxLength(m_index+off); }

  // NOTE: GetMinLoc() could return -1 if you're at m_loc=-4 of a gap(3,5), eg.
  int GetMinLoc() const
    { return ((m_loc >= 0) ? m_loc : Minimum() + m_loc); }
  int GetMaxLoc() const
    { return ((m_loc >= 0) ? m_loc : Maximum() + m_loc); }

  // No error checking done on these values.  For example, SetKmer
  // will happily set m_loc to a range outside what m_index allows.
  void SetPath(const KmerPath& path) { mp_path = &path; }
  void SetIndex(int index) { m_index=index; }
  void SetLoc(int loc) {m_loc = loc; }
  void SetKmer(kmer_id_t k) { m_loc = k-Start(); }


  // These return the # of kmers skipped, or garbage for a gap.
  // (There could be Min and Max versions of them; implement if needed.)
  int IncrementInterval( ) {
    if( atLast() ) return 0;
    int d = mp_path->Length(m_index) - m_loc + 1;
    m_index++;
    m_loc=0;
    return d;
  }
  int DecrementInterval( ) {
    if( atFirst() ) return 0;
    int d = m_loc + 1;
    m_index--;
    if( isSeq() )
      m_loc = mp_path->Length(m_index)-1;
    else
      m_loc = -1; // the right end of a gap
    return d;
  }

  Bool IncrementHaltAtGap( int i );  // returns False if it hits a gap or end
  Bool IncrementMinGap( int i );     // returns False if it falls off the end
  Bool IncrementMaxGap( int i );     // returns False if it falls off the end

  friend inline KmerPathLoc operator+( const KmerPathLoc& loc, int i )
  {    KmerPathLoc answer = loc;
       ForceAssert( answer.IncrementHaltAtGap(i) );
       return answer;    }

  KmerPathLoc& operator+=( nkmers_t i )
  {
    IncrementHaltAtGap( i );
    return *this;
  }

  KmerPathLoc& operator-=( nkmers_t i )
  {
    IncrementHaltAtGap( -i );
    return *this;
  }

  friend inline KmerPathLoc operator-( const KmerPathLoc& loc, int i )
  {    KmerPathLoc answer = loc;
       ForceAssert( answer.IncrementHaltAtGap(-i) );
       return answer;    }

  // Step one to the left/right.
  // Unlike the above, these will step into a gap(0,0), not over it
  bool Decrement()
  { return ((isGap() || atIntervalLeft())
	    ? DecrementInterval()
	    : IncrementHaltAtGap(-1)); }
  bool Increment()
  { return ((isGap() || atIntervalRight())
	    ? IncrementInterval()
	    : IncrementHaltAtGap(+1)); }



  Bool FindKmerBeforeGap( kmer_id_t gapmin, kmer_id_t gapmax,
			  kmer_id_t kmer, vec<KmerPathLoc>& ans );
  Bool FindKmerAfterGap( kmer_id_t gapmin, kmer_id_t gapmax,
			 kmer_id_t kmer, vec<KmerPathLoc>& ans );

  friend ostream& operator<<(ostream& out, KmerPathLoc loc) {
    if( loc.isSeq() )
      out << "kmer " << loc.GetKmer();
    else
      out << "m_loc " << loc.GetLoc();
    out << " in seg " << loc.GetIndex()	<< " of " << loc.GetPath();
    return out;
  }


  friend bool operator== ( const KmerPathLoc& loc1, const KmerPathLoc& loc2 ) {
    return( loc1.mp_path == loc2.mp_path
	    && loc1.m_index == loc2.m_index
	    && loc1.m_loc == loc2.m_loc );
  }

  friend bool operator!= ( const KmerPathLoc& loc1, const KmerPathLoc& loc2 )
  { return( !(loc1 == loc2) ); }

  friend bool operator<( const KmerPathLoc& loc1, const KmerPathLoc& loc2 ) {
    return( loc1.mp_path < loc2.mp_path ||
	    ( loc1.mp_path == loc2.mp_path &&
	      ( loc1.m_index < loc2.m_index ||
		( loc1.m_index == loc2.m_index &&
		  ( loc1.m_loc < loc2.m_loc ) ) ) ) );
  }

  friend bool operator<=( const KmerPathLoc& loc1, const KmerPathLoc& loc2 ) {
    return loc1 < loc2  ||  loc1 == loc2;
  }

  friend bool operator>( const KmerPathLoc& loc1, const KmerPathLoc& loc2 ) {
    return ! ( loc1 <= loc2 );
  }

  friend bool operator>=( const KmerPathLoc& loc1, const KmerPathLoc& loc2 ) {
    return ! ( loc1 < loc2 );
  }

 private:
  // points at RPI mp_path[m_index], k-mer number m_loc
  const KmerPath* mp_path;
  int m_index;
  int m_loc;

};  // class KmerPathLoc



// Global functions relating to KmerPathLoc:
pair<KmerPathLoc,KmerPathLoc>
CreateAlignedLocs( const KmerPath& p1, const KmerPath& p2,
		   int ind1, int ind2 );
// Number of kmers apart; returns 0 if loc1/2 not in the same gap-free path seg
int operator-(KmerPathLoc loc1, KmerPathLoc loc2);
// Distance, crossing gaps using min or max length.
// (NB: if locA is to the right of locB, then DistMax <= DistMin < 0 )
int DistMin(KmerPathLoc locA, KmerPathLoc locB);
int DistMax(KmerPathLoc locA, KmerPathLoc locB);
// CAREFUL: you'll need to add 1 to get a kmer count for a closed interval!

// Count kmers only; gaps count as zero.  This acts on closed intervals,
// ie both endpoints (if they are on kmers, not in gaps!) count.
int KmersInInterval(KmerPathLoc locA, KmerPathLoc locB);


// These push the locs as far as possible over perfectly
// matching stretches of k-mers (not gaps!).  Returns false
// if it had to stop due to k-mer mismatch, ie improper alignment,
// true if it stops due to a gap or beginning/end of read.
bool ScanLeftPerfectMatch(KmerPathLoc& loc1, KmerPathLoc& loc2);
bool ScanRightPerfectMatch(KmerPathLoc& loc1, KmerPathLoc& loc2);
// Same, but they hop identical gaps as well.
// Return true only if stop was due to beginning/end of one sequence.
// These ignore the m_loc of a KmerPathLoc pointing into a gap.
bool ScanLeftPerfectMatchGaps(KmerPathLoc& loc1, KmerPathLoc& loc2);
bool ScanRightPerfectMatchGaps(KmerPathLoc& loc1, KmerPathLoc& loc2);


// Returns true if there is a perfect match from the given locations
// to the end of at least one of the two paths in the stated
// direction, and moves the locations to the last point where the
// paths match.  Returns false if not, leaving the locations in an
// indeterminate state.
// The segments pointed to by the arguments must overlap (else they assert),
// but they ignore the specific kmer pointed to within the segment.
bool IsPerfectMatchToLeft( KmerPathLoc& thisLoc,
                           KmerPathLoc& otherLoc );
bool IsPerfectMatchToRight( KmerPathLoc& thisLoc,
                            KmerPathLoc& otherLoc );
// Returns true if there is a perfect end-to-end match starting
// from the given pair of locations.
bool IsPerfectMatch( KmerPathLoc thisLoc, KmerPathLoc otherLoc );

// Returns true if the given paths have any perfect end-to-end match
bool HavePerfectMatch( const KmerPath& path1, const KmerPath& path2 );


// These methods of KmerPath can't appear until after the full
// definition of KmerPathLoc, so here they are:
inline KmerPathLoc KmerPath::Begin() const
  { return KmerPathLoc(*this,0,0); }
inline KmerPathLoc KmerPath::End() const
  { return KmerPathLoc(*this, size()-1, back().Length()-1 ); }
inline void KmerPath::CopyHead( const KmerPathLoc& loc, KmerPath& ans ) const
  { CopySubpath( Begin(), loc, ans ); }
inline void KmerPath::CopyTail( const KmerPathLoc& loc, KmerPath& ans ) const
  { CopySubpath( loc, End(), ans ); }
inline void KmerPath::CopyHeadNoLastKmer(
     const KmerPathLoc& loc, KmerPath& ans ) const
  { CopySubpathNoLastKmer( Begin(), loc, ans ); }
inline void KmerPath::CopyTailNoFirstKmer(
     const KmerPathLoc& loc, KmerPath& ans ) const
  { CopySubpathNoFirstKmer( loc, End(), ans ); }


// FuncDecl: ProperOverlapExt
//
// Given a fixed position at which two kmer paths without gaps
// agree, determine if it extends to a proper overlap between them.  The fixed
// positions are defined by the indices of the segments in the kmer paths, and by
// the common position, except that the answer can be computed without knowing this
// common position.
//
// Don't call this with KmerPaths having gaps.
Bool ProperOverlapExt( const KmerPath& p1, const KmerPath& p2, int ind1, int ind2 );

// CreateDatabase: Input paths and (maybe) paths_rc.  Create a vec<TAG> that
// consists of every KmerPathInterval in these paths, sorted by kmer_id, but
// tagged with the interval's original location in the paths.  This pathsdb can
// then be used for quick lookups of a kmer_id via Contains(...).
template<class TAG> void CreateDatabase(
     const vecKmerPath& paths, vec<TAG>& pathsdb );
template<class TAG> void CreateDatabase(
     const vecKmerPath& paths, const vecKmerPath& paths_rc, vec<TAG>& pathsdb );

// Find which reads contain kmers NOT in the reference read path database.
// WARNING: Make surethe reads and the reference database use the same kmer
// numbering scheme!  If this function produces an unexpectedly high number of
// novel reads, the kmers probably don't match up.
template <class TAG>
void
MarkReadsWithNovelKmers( const vecKmerPath & reads, const vec<TAG>& reference,
			 vec<Bool>& answer );

// As above, but uses adjacency graph for more stringent testing
template <class TAG>
void
MarkReadsWithNovelKmers2( const vecKmerPath & reads, const vec<TAG>& reference,
			  const digraph& adj_graph, vec<Bool>& answer );


// FuncDecl: SubContig
//
// Given a KmerPath p, test whether it's genomic as follows:
// determine if it can be built up from a contig
// consisting of paths (in paths, paths_rc) referenced by xpathsdb.  Note that
// xpathsdb may have been built from only some of the KmerPaths in paths and
// paths_rc.  Ignore any paths in xpathsdb with ReadId()s in pPathsToIgnore, which
// should be sorted, and (if specified) the path with ReadId() equal to pathToIgnore.
Bool SubContig( const KmerPath& p, const vecKmerPath& paths,
                const vecKmerPath& paths_rc, const vec<tagged_rpint>& xpathsdb,
                const HashSimple * pPathsToIgnore = 0, read_id_t pathToIgnore = NULL_READ_ID );

/**
   Class: vecKmerPathIndex

   An index of a <vecKmerPath>, that lets you quickly go from position along a
   path to the segment containing the position.  Normally used to take
   a vector of <genome paths> representing the genome parts, and quickly
   go from genome part id and position within that part to the kmer.
*/
class vecKmerPathIndex {
 public:

  vecKmerPathIndex( const vecKmerPath& _gpaths, const vecKmerPath& _gpaths_rc );


  // Method: FindLoc
  //
  // Given a position in the path, find the segment and the location within the segment
  // for the kmer starting at that position.
  KmerPathLoc FindLoc( longlong pathId, int posInPath, Bool orient ) const;

 private:
  const vecKmerPath& gpaths, gpaths_rc;

  // Field: segmentStarts
  // For each path, for each segment, the position (base number in the path) at which
  // that segment's first kmer starts.
  // For genome path of id i, segmentStarts[i] is an array whose j'th element gives the
  // starting position (in bases) of the j'th segment of the path.
  vec< vec< genome_part_pos_t > > segmentStarts, segmentStarts_rc;

};  // class KmerPathIndex

// Semantic type: path_id_t
// Identifies a <path> by its index in an array of paths.
SemanticType( basevec_id_t, path_id_t );


// Logical type: unipath_id_t
// Identifies a <unipath> by its index in an array of unipaths.
SemanticType( path_id_t, unipath_id_t );

// Logical type: unipath_interval_id_t
// Index of a <unipath> interval in a <path interval database>
// for a set of unipaths.
SemanticType( path_interval_id_t, unipath_interval_id_t );


#endif

// Synonyms: Various synonyms
//    kmer path - See <KmerPath>
//    read path - See <KmerPath>
//    segment - See <KmerPathInterval>

// Term: path
//
// When used without qualifiers, "path" refers to "read path".



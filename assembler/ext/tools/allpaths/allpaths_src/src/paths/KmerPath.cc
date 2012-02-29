///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "system/MemTracker.h"

// Build a KmerPath from the String of its <<.
// Halts silently at any syntax errors.
KmerPath::KmerPath( String s ) {
  s.GlobalReplaceBy(" ", "");
  if(s[0]=='{') 
    s=s.After("{");
  while(!s.empty()) {
    if( s[0]=='[' ) {
      this->AddSegment( s.Between("[","-").Int(), s.Between("-","]").Int() );
      s = s.After("]");
    }
    else if( s[0]=='(' ) {
      this->AddGap( s.Between("gap","-").Int(), s.Between("-",")").Int() );
      s = s.After(")");
    }
    else break;
  }
}

// Build a KmerPath from a vec of kmer ids
// It does NOT check to see if the path is valid (yet)
KmerPath::KmerPath(const vec<kmer_id_t>& path ) {
  for (size_t i = 0; i < path.size(); ++i) 
    this->AddSegment( path[i], path[i] );
  this->Canonicalize();
}

/**
 * BasePosToSegmentPos
 *
 * Input: base position in the bvec encoded by this KmerPath. Output:
 * corresponding KmerPathLoc of kmer starting at base_pos.
 *
 * Set index and loc to -1 on error.
 */
KmerPathLoc KmerPath::BasePosToSegmentPos( int base_pos ) const
{
  int curr_len = 0;
  for (int ii=0; ii<this->NSegments( ); ii++) {
    longlong start_kpi = this->Segment( ii ).Start( );
    longlong stop_kpi =  this->Segment( ii ).Stop( );
    longlong quasi_len = stop_kpi - start_kpi;
    if ( curr_len <= base_pos && base_pos <= curr_len + quasi_len )
      return KmerPathLoc( *this, ii, base_pos - curr_len );
    curr_len += 1 + quasi_len;
  }
  
  // Not found.
  return KmerPathLoc( *this, -1, -1 );
}

// Methods and friend functions of class KmerPath:

void KmerPath::Canonicalize( ) {
  int currSeg = 0;
  int nextSeg = 1;
  int numSegs = this->size();
  while ( nextSeg < numSegs )
  {
    if ( this->isGap(currSeg) && this->isGap(nextSeg) )
    {
      // Join consecutive gaps.
      this->SetMinimum(currSeg, this->Minimum(currSeg) + this->Minimum(nextSeg) );
      this->SetMaximum(currSeg, this->Maximum(currSeg) + this->Maximum(nextSeg) );
    }
    else if ( this->isSeq(currSeg) && this->isSeq(nextSeg) &&
              this->Stop(currSeg) + 1 == this->Start(nextSeg) )
    {
      // Join consecutive mergable sequence segments.
      this->SetStop(currSeg, this->Stop(nextSeg) );
    }
    else
    {
      // Copy unmergeable segment to just after current segment, and
      // make that the current segment.
      ++currSeg;
      (*this)[currSeg] = (*this)[nextSeg];
    }
    
    // Check the following segment.
    ++nextSeg;
  }

  // Remove superfluous segments.
  this->resize( currSeg+1 );
}

void KmerPath::AddSegment( KmerPathInterval rpi ) {
  if (size()==0) // then don't ask about back()!
    push_back( rpi );
  // Concatenate abutting intervals of sequence:
  else if (rpi.isSeq() && back().isSeq() && rpi.Start() == 1+back().Stop())
    back().Set( back().Start(), rpi.Stop(), rpi.isGap() );
  // Combine consecutive gaps (maybe never used):
  else if (rpi.isGap() && back().isGap())
    back().Set( back().Start() + rpi.Start(), 
		back().Stop() + rpi.Stop(), rpi.isGap() );
  else
    push_back( rpi );
}

void KmerPath::AddSegmentNoConcatenate( KmerPathInterval rpi ) {
    push_back( rpi );
}

// Defaults of begin=0, end=-1
void KmerPath::Append( const KmerPath& other, int begin, int end ) {
  if (end<0) end += other.NSegments();
  this->reserve( this->size() + end-begin+1 + 1 );
  for(int i=begin; i<=end; i++)
    this->AddSegment( other.Segment(i) );
}


// Defaults of begin=0, end=-1
// Don't call this on a KmerPath starting with a gap:
void KmerPath::AppendNoFirstKmer( const KmerPath& other, int begin, int end ) {
  Assert( other.isSeq(begin) );
  if (end<0) end += other.NSegments();
  this->reserve( this->size() + end-begin+1 + 1 );
  if( other.Length(begin) > 1 )
    this->AddSegment( other.Start(begin)+1, other.Stop(begin),
		      other.isGap(begin) );
  Append(other, begin+1, end);
}



// The answer is appended to ans, so you can build up a path
// with successive calls to Subpath and SubpathNoFirstKmer:
void KmerPath::CopySubpath( KmerPathLoc loc1, KmerPathLoc loc2, 
			    KmerPath& ans ) const {
  // loc1 and loc2 must point to *this
  ForceAssert( loc1.GetPathPtr() == this );
  ForceAssert( loc2.GetPathPtr() == this );

  int i1=loc1.GetIndex(), i2=loc2.GetIndex();
  // If they are (even possibly) in the wrong order, return (the empty subpath):
  if( i1>i2 || ( i1==i2 && loc1.GetMaxLoc() > loc2.GetMinLoc()) )
    return;

  // Handle the same interval case:
  if( i1==i2 ) {
    if( loc1.isSeq() )
      ans.AddSegment( loc1.GetKmer(), loc2.GetKmer() );
    else
      ans.AddGap( loc2.GetMinLoc()-loc1.GetMaxLoc() + 1,
		  loc2.GetMaxLoc()-loc1.GetMinLoc() + 1 );
    return;
  }

  ans.reserve( ans.size() + i2 - i1 + 1 );

  // The general case:
  // Copy the first KmerPathInterval:
  if( loc1.isSeq() )
    ans.AddSegment( loc1.GetKmer(), loc1.Stop() );
  else  // This yeilds a gap with no wiggle room if m_loc measuring from right
    ans.AddGap( max((longlong)0, loc1.Minimum()-loc1.GetMinLoc()),
		loc1.Maximum()-loc1.GetMaxLoc() );
  // Copy the middle ones:
  for ( int i = i1+1; i < i2; ++i )
    ans.push_back( this->Segment(i) );
  // Copy the last one:
  if( loc2.isSeq() )
    ans.AddSegment( loc2.Start(), loc2.GetKmer() );
  else  // This yeilds a gap with no wiggle room if m_loc measuring from left
    ans.AddGap( max(0, loc2.GetMinLoc()), loc2.GetMaxLoc() );
}


void KmerPath::CopySubpathNoFirstKmer( KmerPathLoc loc1, KmerPathLoc loc2, 
				       KmerPath& ans ) const {
  // loc1 and loc2 must point to *this
  ForceAssert( loc1.GetPathPtr() == this );
  ForceAssert( loc2.GetPathPtr() == this );

  // I think the NoFirstKmer version should choke if loc1 is on a gap.
  // (But perhaps this is wrong and it should just shorten the gap by one?)
  ForceAssert( ! loc1.isGap() );

  int i1=loc1.GetIndex(), i2=loc2.GetIndex();
  // If they are in the wrong order, return (the empty subpath):
  if( i1>i2 || ( i1==i2 && loc1.GetLoc() >= loc2.GetLoc()) ) // NoFirstKmer >=
    return;

  // Handle the same interval case:
  if( i1==i2 ) {
    ans.AddSegment( loc1.GetKmer()+1, loc2.GetKmer() );
    return;
  }

  ans.reserve( ans.size() + i2 - i1 + 1 );

  // The general case:
  // Copy the first KmerPathInterval:
  if( ! loc1.atIntervalRight() )
    ans.AddSegment( loc1.GetKmer()+1, loc1.Stop() );
  // Copy the middle ones:
  if ( i1+1 < i2 )
    ans.AddSegment( this->Segment(++i1) );
  for ( int i = i1+1; i < i2; ++i )
    ans.push_back( this->Segment(i) );
  // Copy the last one:
  if( loc2.isSeq() )
    ans.AddSegment( loc2.Start(), loc2.GetKmer() );
  else  // This yeilds a gap with no wiggle room if m_loc measuring from left
    ans.AddGap( max(0, loc2.GetMinLoc()), loc2.GetMaxLoc() );
}

void KmerPath::CopySubpathNoLastKmer( KmerPathLoc loc1, KmerPathLoc loc2, 
				      KmerPath& ans ) const {
  // I think the NoLastKmer version should choke if loc2 is on a gap.
  // (But perhaps this is wrong and it should just shorten the gap by one?)
  ForceAssert( ! loc2.isGap() );

  CopySubpath( loc1, loc2, ans );
  
  // pop the last kmer:
  if( loc2.atIntervalLeft() || loc1 == loc2 )
    // ans.pop_back() -- oops, can't do that:
    ans.SetNSegments( ans.NSegments()-1 );
  else
    ans.SetStop( ans.NSegments()-1, ans.Stop(ans.NSegments()-1)-1 );

}

// Here's the hard one: copy the path [loc1,loc2] into the space
// of a gap with size given_min to given_max.  In practice, all this
// means is that if we copy a gap which couldn't possibly attain its
// maximum length in its new home -- even with all other gaps at their
// minima -- then reduce its max length appropriately (&likewise for min).

void KmerPath::CopySubpathAdjustGaps( KmerPathLoc loc1, KmerPathLoc loc2, 
				      int given_min, int given_max, 
				      KmerPath& ans ) const {

  // loc1 and loc2 must point to *this
  ForceAssert( loc1.GetPathPtr() == this );
  ForceAssert( loc2.GetPathPtr() == this );

  int i1=loc1.GetIndex(), i2=loc2.GetIndex();
  // If they are in the wrong order, return (the empty subpath):
  if( i1>i2 || ( i1==i2 && loc1.GetMaxLoc() > loc2.GetMinLoc()) )
    return;

  // Make sure the interval fits:
  int orig_min = DistMin(loc1,loc2)+1, orig_max = DistMax(loc1,loc2)+1;

  // Remember those +1's for closed intervals!
  ForceAssertLe(orig_min, given_max);
  ForceAssertLe(given_min, orig_max);
  // So these measures of wiggle room must be >=0:
  int shrink_room = orig_max - given_min;
  int grow_room = given_max - orig_min;

  // Handle the same interval case:
  if( i1==i2 ) {
    if( loc1.isSeq() )
      ans.AddSegment( loc1.GetKmer(), loc2.GetKmer() );
    else
      ans.AddGap( max( given_min, loc2.GetMinLoc()-loc1.GetMaxLoc() + 1 ),
		  min( given_max, loc2.GetMaxLoc()-loc1.GetMinLoc() + 1 ) );
    return;
  }

  // The general case:

  ans.reserve( ans.size() + i2-i1+1 + 1 );

  // Copy the first KmerPathInterval
  if( loc1.isSeq() )
    ans.AddSegment( loc1.GetKmer(), loc1.Stop() );
  else {
    int mymin = max((longlong)0, loc1.Minimum()-loc1.GetMinLoc());
    int mymax = loc1.Maximum()-loc1.GetMaxLoc();
    ans.AddGap( max( mymin, mymax - shrink_room ),
		min( mymax, mymin + grow_room ) );
  }
  // Copy the middle ones.  Adjust gaps as needed.
  for(int i=i1+1; i <= i2-1; i++) {
    if( isSeq(i) )
      ans.AddSegment( Segment(i) );
    else
      ans.AddGap( max( Minimum(i), Maximum(i)-shrink_room ),
		  min( Maximum(i), Minimum(i)+grow_room ) );
  }
  // Copy the last one 
  if( loc2.isSeq() )
    ans.AddSegment( loc2.Start(), loc2.GetKmer() );
  else {
    int mymin = max(0, loc2.GetMinLoc() + 1);   // +1's for closed intervals:
    int mymax = loc2.GetMaxLoc() + 1;
    ans.AddGap( max( mymin, mymax - shrink_room ),
		min( mymax, mymin + grow_room ) );
  }
}


Bool IsSubpathAnchoredLeft( const KmerPath& p1, const KmerPath& p2 )
{    int n1 = p1.NSegments( ), n2 = p2.NSegments( );
     if ( n1 > n2 ) return False;
     if ( n1 == 0 ) return True;
     for ( int i = n1 - 2; i >= 0; i-- )
          if ( p1.Segment(i) != p2.Segment(i) ) return False;
     KmerPathInterval s1 = p1.Segment( n1 - 1 ), s2 = p2.Segment( n1 - 1 );
     if ( s1.isGap( ) )
     {    if ( !s2.isGap( ) ) return False;
          return s1.Maximum( ) <= s2.Maximum( );    }
     if ( s2.isGap( ) ) return False;
     if ( s1.Start( ) != s2.Start( ) ) return False;
     return s1.Stop( ) <= s2.Stop( );    }

void PrintFolded( ostream& out, const KmerPath& p )
{   out << "{";
    for( int i = 0; i < p.NSegments( ); i++ )
    {    if ( i > 0 && i % 3 == 0 ) cout << "\n";
         out << p.Segment(i);    }
    out << "}\n";    }

ostream& operator<<(ostream& out, const KmerPath& p)
{   out << "{";
    for( int i = 0; i < p.NSegments( ); i++ )
         out << p.Segment(i);
    return out << "}";    }

void KmerPath::AppendToDatabase( vec<tagged_rpint>& segs, int i ) const
{    for ( int j = 0; j < NSegments( ); j++ )
          Segment(j).AppendToDatabase( segs, i, j );    }

void KmerPath::AppendToDatabase( vec<big_tagged_rpint>& segs, int i ) const
{    for ( int j = 0; j < NSegments( ); j++ )
          Segment(j).AppendToDatabase( segs, i, j );    }

void KmerPath::AppendToDatabase( vec<new_tagged_rpint>& segs, int i ) const
{    for ( int j = 0; j < NSegments( ); j++ )
          Segment(j).AppendToDatabase( segs, i, j );    }

void KmerPath::Reverse( )
{    KmerPath r;
     for ( int j = NSegments( ) - 1; j >= 0; j-- )
     {    const KmerPathInterval& I = Segment(j);
          if ( I.isGap( ) ) r.AddSegment(I);
          else if ( is_palindrome( I.Start( ) ) )
          {    for ( longlong u = I.Stop( ); u >= I.Start( ); u-- )
                    r.AddSegment( u, u );    }
          else r.AddSegment( flip_kmer( I.Stop( ) ), flip_kmer( I.Start( ) ) );    }
     *this = r;    }

Bool KmerPath::Proper( ) const
{    if ( NSegments( ) == 0 ) return False; 
     if ( Segment(0).isGap( ) || LastSegment( ).isGap( ) ) return False;
     for ( int i = 1; i < NSegments( ); i++ )
     {    if ( Segment(i).isSeq( ) && Segment(i-1).isSeq( ) 
               && Segment(i).Start( ) == 1 + Segment(i-1).Stop( ) )
          {    return False;    }
          if ( Segment(i).isGap( ) && Segment(i-1).isGap( ) ) return False;    }
     return True;    }

void BinaryWrite( int fd, const KmerPath& p )
{    int n = p.NSegments( );
     WriteBytes( fd, &n, sizeof(int) );
     WriteBytes( fd, &p.Segment(0), n * sizeof(KmerPathInterval) );    }

void BinaryRead( int fd, KmerPath& p )
{    int n;
     ReadBytes( fd, &n, sizeof(int) );
     p.SetNSegments(n);
     ReadBytes( fd, &p.Segment(0), n * sizeof(KmerPathInterval) );    }
  


// Methods of class KmerPathLoc:

// Method: IncrementHaltAtGap
//
// Returns False if it hits a gap or the end of the path.
// You can tell the difference by checking atFirst() or atLast().
Bool KmerPathLoc::IncrementHaltAtGap( int i ) {
  if( m_index<0 || m_index>=mp_path->NSegments() || isGap() )
    return False;
  m_loc += i;
  for(;;) {
    if( m_loc >= Length() ) {
      m_index++;
      if( m_index >= mp_path->NSegments() || isGap() ) {
	// oops; undo:
	m_index--;
	m_loc = Length() - 1;
	return False;
      } else
	m_loc -= Length(-1);
    } else if( m_loc < 0 ) {
      m_index--;
      if( m_index < 0 || isGap() ) {
	// oops; undo:
	m_index++;
	m_loc = 0;
	return False;
      } else
	m_loc += Length();
    } else
      return True;
  }
}

// Increment using Min gap sizes.
Bool KmerPathLoc::IncrementMinGap( int i ) {
  if( m_index<0 || m_index>=mp_path->NSegments() )
    return False;
  Bool from_right = (m_loc<0);
  if( from_right ) m_loc += MinLength();
  m_loc += i;
  for(;;) {
    if( m_loc >= MinLength() ) {
      from_right = False;
      m_index++;
      if( m_index >= mp_path->NSegments() ) {
	// oops; undo:
	m_index--;
	m_loc = MinLength() - 1;
	return False;
      } else
	m_loc -= MinLength(-1);
    } else if( m_loc < 0 ) {
      from_right = True;
      m_index--;
      if( m_index < 0 ) {
	// oops; undo:
	m_index++;
	m_loc = 0; // a KmerPath should never start with a gap
	return False;
      } else
	m_loc += MinLength();
    } else {
      if( isGap() && from_right ) // set m_loc to -(# gap kmers used)-1
	m_loc -= MinLength();
      return True;
    }
  }
}

// Increment using Max gap sizes.
Bool KmerPathLoc::IncrementMaxGap( int i ) {
  if( m_index<0 || m_index>=mp_path->NSegments() )
    return False;
  Bool from_right = (m_loc<0);
  if( from_right ) m_loc += MaxLength();
  m_loc += i;
  for(;;) {
    if( m_loc >= MaxLength() ) {
      from_right = False;
      m_index++;
      if( m_index >= mp_path->NSegments() ) {
	// oops; undo:
	m_index--;
	m_loc = MaxLength() - 1;
	return False;
      } else
	m_loc -= MaxLength(-1);
    } else if( m_loc < 0 ) {
      from_right = True;
      m_index--;
      if( m_index < 0 ) {
	// oops; undo:
	m_index++;
	m_loc = 0; // a KmerPath should never start with a gap
	return False;
      } else
	m_loc += MaxLength();
    } else {
      if( isGap() && from_right ) // set m_loc to -(# gap kmers used)-1
	m_loc -= MaxLength();
      return True;
    }
  }
}



// Finds all occurrences of a given k-mer in [min,max], where
//   ---:RPL:---------min----max-----
// Stops looking and returns False there's a gap between RPL and max.

Bool KmerPathLoc::FindKmerAfterGap( longlong gapmin, longlong gapmax,
				    longlong kmer, vec<KmerPathLoc>& ans ) {
  if( IncrementHaltAtGap(gapmin+1) == False )
    return atLast();  // hit a gap or end of sequence.

  // now RPL points at the first k-mer that could be after the gap
  int to_search = gapmax-gapmin;
  while( to_search >= 0 ) {
    if( isGap() ) return False;
    if( GetKmer() <= kmer 
	&& kmer <= GetKmer() + to_search
	&& kmer <= Stop() ) {
      ans.push_back(*this);
      ans.back().IncrementHaltAtGap( kmer - GetKmer() );
    }
    if( atLast() ) return True;
    to_search -= IncrementInterval();
  }
  return True;
}
// Same, but moving left:
Bool KmerPathLoc::FindKmerBeforeGap( longlong gapmin, longlong gapmax,
				     longlong kmer, vec<KmerPathLoc>& ans ) {
  if( IncrementHaltAtGap(-gapmin-1) == False )
    return atFirst();  // hit a gap or end of sequence.

  // now RPL points at the first k-mer that could be after the gap
  int to_search = gapmax-gapmin;
  while( to_search >= 0 ) {
    if( isGap() ) return False;
    if( Start() <= kmer 
	&& GetKmer() - to_search <= kmer
	&& kmer <= GetKmer() ) {
      ans.push_back(*this);
      ans.back().IncrementHaltAtGap( kmer - GetKmer() );
    }
    if( atFirst() ) return True;
    to_search -= DecrementInterval();
  }
  return True;
}



// Global functions that act on KmerPathLocs.

pair<KmerPathLoc,KmerPathLoc> 
CreateAlignedLocs( const KmerPath& p1, const KmerPath& p2, 
		   int ind1, int ind2 ) {
  ForceAssert( p1.isSeq(ind1) );
  ForceAssert( p2.isSeq(ind2) );
  ForceAssert( p1.Segment(ind1).Overlaps(p2.Segment(ind2)) );
  longlong k = max( p1.Start(ind1), p2.Start(ind2) );
  return( make_pair( KmerPathLoc(p1, ind1, k-p1.Start(ind1)),
		     KmerPathLoc(p2, ind2, k-p2.Start(ind2)) ) );
}


// Number of kmers you'd need to increment loc1 by to get loc2 (+/-).
// Returns 0 if can't do it without crossing a gap.
int operator-(KmerPathLoc locB, KmerPathLoc locA) {
  if( locB.GetPathPtr() != locA.GetPathPtr() ) return 0;
  if( locB.GetIndex() == locA.GetIndex() )
    return locB.GetLoc() - locA.GetLoc();
  Bool exch = ( locB.GetIndex() < locA.GetIndex() );
  if(exch) swap(locB,locA);
  int d = ( locA.Length()-locA.GetLoc() ) + locB.GetLoc();
  for(int i=1; i<locB.GetIndex()-locA.GetIndex(); i++) {
    if( locA.isGap(i) ) return 0;
    d += locA.Length(i);
  }
  return( exch ? -d : d );
}

int DistMin(KmerPathLoc locA, KmerPathLoc locB) {
  if( locB.GetPathPtr() != locA.GetPathPtr() ) return 0;
  if(locA.isGap() && locA.GetLoc()<0) locA.SetLoc(locA.GetLoc()+locA.Minimum());
  if(locB.isGap() && locB.GetLoc()<0) locB.SetLoc(locB.GetLoc()+locB.Minimum());
  if( locB.GetIndex() == locA.GetIndex() )
    return locB.GetLoc() - locA.GetLoc();
  Bool exch = ( locB.GetIndex() < locA.GetIndex() );
  if(exch) swap(locB,locA);
  int d = ( locA.MinLength()-locA.GetLoc() ) + locB.GetLoc();
  for(int i=1; i<locB.GetIndex()-locA.GetIndex(); i++) {
    if( locA.isGap(i) )
      d += locA.Minimum(i);
    else
      d += locA.Length(i);
  }
  return( exch ? -d : d );
}

int DistMax(KmerPathLoc locA, KmerPathLoc locB) {
  if( locB.GetPathPtr() != locA.GetPathPtr() ) return 0;
  if(locA.isGap() && locA.GetLoc()<0) locA.SetLoc(locA.GetLoc()+locA.Maximum());
  if(locB.isGap() && locB.GetLoc()<0) locB.SetLoc(locB.GetLoc()+locB.Maximum());
  if( locB.GetIndex() == locA.GetIndex() )
    return locB.GetLoc() - locA.GetLoc();
  Bool exch = ( locB.GetIndex() < locA.GetIndex() );
  if(exch) swap(locB,locA);
  int d = ( locA.MaxLength()-locA.GetLoc() ) + locB.GetLoc();
  for(int i=1; i<locB.GetIndex()-locA.GetIndex(); i++) {
    if( locA.isGap(i) )
      d += locA.Maximum(i);
    else
      d += locA.Length(i);
  }
  return( exch ? -d : d );
}

int KmersInInterval(KmerPathLoc locA, KmerPathLoc locB) {
  if( locB.GetPathPtr() != locA.GetPathPtr() ) return 0;
  if( locB.GetIndex() == locA.GetIndex() )
    return (locA.isSeq() ? locB.GetLoc() - locA.GetLoc() + 1 : 0);
  bool exch = ( locB.GetIndex() < locA.GetIndex() );
  if(exch) swap(locB,locA);
  int d=0;
  if( locA.isSeq() ) d += locA.Length() - locA.GetLoc();
  if( locB.isSeq() ) d += locB.GetLoc() + 1;
  for(int i=1; i<locB.GetIndex()-locA.GetIndex(); i++)
    if( locA.isSeq(+i) )
      d += locA.Length(+i);
  return( exch ? -d : d );
}
   
    



bool ScanLeftPerfectMatch(KmerPathLoc& loc1, KmerPathLoc& loc2) {
  if( loc1.GetKmer() != loc2.GetKmer() ) return false;
  while( loc1.Start()==loc2.Start()
	 && loc1.isSeq(-1) && loc2.isSeq(-1)
	 && loc1.Stop(-1)==loc2.Stop(-1) ) {
    loc1.DecrementInterval();
    loc2.DecrementInterval();
  }
  longlong k = max( loc1.Start(), loc2.Start() );
  loc1.SetKmer(k);
  loc2.SetKmer(k);

  // Go through the reasons we might have stopped:

  if( loc1.Start()<loc2.Start() ) return( !loc2.isSeq(-1) );
  if( loc1.Start()>loc2.Start() ) return( !loc1.isSeq(-1) );

  if( !loc1.isSeq(-1) || !loc2.isSeq(-1) ) return( true );

  if( loc1.Stop(-1) != loc2.Stop(-1) ) return( false );

  // One of those must be true; we should never get here:
  ForceAssert(0==1);
  return(false);
}

bool ScanRightPerfectMatch(KmerPathLoc& loc1, KmerPathLoc& loc2) {
  if( loc1.GetKmer() != loc2.GetKmer() ) return false;
  while( loc1.Stop()==loc2.Stop()
	 && loc1.isSeq(+1) && loc2.isSeq(+1)
	 && loc1.Start(+1)==loc2.Start(+1) ) {
    loc1.IncrementInterval();
    loc2.IncrementInterval();
  }
  longlong k = min( loc1.Stop(), loc2.Stop() );
  loc1.SetKmer(k);
  loc2.SetKmer(k);

  // Go through the reasons we might have stopped:

  if( loc1.Stop()<loc2.Stop() ) return( !loc1.isSeq(+1) );
  if( loc1.Stop()>loc2.Stop() ) return( !loc2.isSeq(+1) );

  if( !loc1.isSeq(+1) || !loc2.isSeq(+1) ) return( true );

  if( loc1.Start(+1) != loc2.Start(+1) ) return( false );

  // One of those must be true; we should never get here:
  ForceAssert(0==1);
  return(false);
}

bool ScanRightPerfectMatchGaps(KmerPathLoc& loc1, KmerPathLoc& loc2) {
  // make sure we start with a match:
  if( ! (loc1.isSeq() && loc2.isSeq() && loc1.GetKmer() == loc2.GetKmer())
      &&
      ! (loc1.isGap() && loc2.isGap() && loc1.GetSegment() == loc2.GetSegment()) )
    return false;

  while( // The current ones are the same (on the right, if seq)
        ((loc1.isGap() && loc2.isGap() && loc1.GetSegment() == loc2.GetSegment()) ||
         (loc1.isSeq() && loc2.isSeq() && loc1.Stop() == loc2.Stop()))
	&&
	// The intervals to the right match (on the left, if seq)
	((loc1.isGap(+1) && loc2.isGap(+1) && loc1.GetSegment(+1) == loc2.GetSegment(+1)) ||
	 (loc1.isSeq(+1) && loc2.isSeq(+1) && loc1.Start(+1)==loc2.Start(+1)))
	) {
    loc1.SetIndex(loc1.GetIndex()+1);
    loc2.SetIndex(loc2.GetIndex()+1);
  }
  // Set locs to last matching kmer
  if( loc1.isSeq() && loc2.isSeq() ) {
    longlong k = min( loc1.Stop(), loc2.Stop() );
    loc1.SetKmer(k);
    loc2.SetKmer(k);
  }

  // Return true only if we stopped due to running out of sequence
  return( loc1.atEnd() || loc2.atEnd() );
}

bool ScanLeftPerfectMatchGaps(KmerPathLoc& loc1, KmerPathLoc& loc2) {
  // make sure we start with a match:
  if( ! (loc1.isSeq() && loc2.isSeq() && loc1.GetKmer() == loc2.GetKmer())
      &&
      ! (loc1.isGap() && loc2.isGap() && loc1.GetSegment() == loc2.GetSegment()) )
    return false;

  while( // The current ones are the same (on the left, if seq)
        ((loc1.isGap() && loc2.isGap() && loc1.GetSegment() == loc2.GetSegment()) ||
         (loc1.isSeq() && loc2.isSeq() && loc1.Start() == loc2.Start()))
	&&
	// The intervals to the left match (on the right, if seq)
	((loc1.isGap(-1) && loc2.isGap(-1) && loc1.GetSegment(-1) == loc2.GetSegment(-1)) ||
	 (loc1.isSeq(-1) && loc2.isSeq(-1) && loc1.Stop(-1)==loc2.Stop(-1)))
	) {
    loc1.SetIndex(loc1.GetIndex()-1);
    loc2.SetIndex(loc2.GetIndex()-1);
  }
  // Set locs to last matching kmer
  if( loc1.isSeq() && loc2.isSeq() ) {
    longlong k = max( loc1.Start(), loc2.Start() );
    loc1.SetKmer(k);
    loc2.SetKmer(k);
  }

  // Return true only if we stopped due to running out of sequence
  return( loc1.atBegin() || loc2.atBegin() );
}


// Returns true if there is a perfect match to the left of the given
// locations, and moves the locations to the last point where the
// paths match.  Returns false if not, leaving the locations in an
// indeterminate state.
// The segments pointed to by the arguments must overlap (else it asserts),
// but it ignores the specific kmer pointed to within the segment.
bool IsPerfectMatchToLeft( KmerPathLoc& thisLoc,
                           KmerPathLoc& otherLoc )
{

  ForceAssert( thisLoc.GetSegment().Overlaps( otherLoc.GetSegment() ) );

  const KmerPath* p_thisPath = thisLoc.GetPathPtr();
  int thisSeg = thisLoc.GetIndex();
  
  const KmerPath* p_otherPath = otherLoc.GetPathPtr();
  int otherSeg = otherLoc.GetIndex();

  // If we're not at the beginning of either path...
  if ( otherSeg > 0 && 
       thisSeg > 0 )
  {  
    // ...and we can't jump to the previous segment, fail.
    if ( p_otherPath->Start( otherSeg ) != p_thisPath->Start( thisSeg ) )
      return false;
    
    // ...but if we can jump back...
    else
    {
      // ...walk back on both paths till you hit a difference or the
      // first segment of one or both paths.
      int minNumSegsToLeft = min( otherSeg, thisSeg );
      
      while ( minNumSegsToLeft > 0 &&
              p_otherPath->Segment( --otherSeg ) == p_thisPath->Segment( --thisSeg ) )
        --minNumSegsToLeft;
      
      // If we hit a difference before the final segment OR we're at
      // the final segment and one or the other path has a gap there
      // or they're both sequence there but they stop at different
      // kmers, then the match isn't perfect.
      
      if ( minNumSegsToLeft > 1 ||
           ( minNumSegsToLeft == 1 &&
             ( p_otherPath->isGap( otherSeg ) ||
               p_thisPath->isGap( thisSeg ) ||
               p_otherPath->Stop( otherSeg ) != p_thisPath->Stop( thisSeg ) ) ) )
        return false;
    }
  }

  longlong kmer = max( p_otherPath->Start( otherSeg ), p_thisPath->Start( thisSeg ) );

  otherLoc.SetIndex( otherSeg );
  otherLoc.SetKmer( kmer );
  thisLoc.SetIndex( thisSeg );
  thisLoc.SetKmer( kmer );

  return ( otherLoc.atBegin() || thisLoc.atBegin() );
}


// Returns true if there is a perfect match to the right of the given
// locations, and moves the locations to the last point where the
// paths match.  Returns false if not, leaving the locations in an
// indeterminate state.
// The segments pointed to by the arguments must overlap (else it asserts),
// but it ignores the specific kmer pointed to within the segment.
bool IsPerfectMatchToRight( KmerPathLoc& thisLoc,
                            KmerPathLoc& otherLoc )
{
  ForceAssert( thisLoc.GetSegment().Overlaps( otherLoc.GetSegment() ) );

  const KmerPath* p_thisPath = thisLoc.GetPathPtr();
  int thisSeg = thisLoc.GetIndex();

  const KmerPath* p_otherPath = otherLoc.GetPathPtr();
  int otherSeg = otherLoc.GetIndex();

  // If we're not at the end of one sequence or the other...
  if ( ! thisLoc.atLast() &&
       ! otherLoc.atLast() )
  {
    // ...and we can't jump to the next segment, fail.
    if ( p_otherPath->Stop( otherSeg ) != p_thisPath->Stop( thisSeg ) )
      return false;
    
    // ...but if we can jump...
    else
    {
      // ...advance on both paths till you hit a difference or the end
      // of one or both paths.
      int minNumSegsToRight = min( p_otherPath->NSegments() - otherSeg,
                                   p_thisPath->NSegments() - thisSeg ) - 1;
      
      while ( minNumSegsToRight > 0 &&
              p_otherPath->Segment( ++otherSeg ) == p_thisPath->Segment( ++thisSeg ) )
        --minNumSegsToRight;
      
      // If we hit a difference before the last segment, or we're at
      // the last segment of one or both paths and they're
      // incompatible, fail.
      if ( minNumSegsToRight > 1 ||
           ( minNumSegsToRight == 1 &&
             ( p_otherPath->isGap( otherSeg ) ||
               p_thisPath->isGap( thisSeg ) ||
               p_otherPath->Start( otherSeg ) != p_thisPath->Start( thisSeg ) ) ) )
        return false;
    }
  }
            
  longlong kmer = min( p_otherPath->Stop( otherSeg ), p_thisPath->Stop( thisSeg ) );
  otherLoc.SetIndex( otherSeg );
  otherLoc.SetKmer( kmer );
  thisLoc.SetIndex( thisSeg );
  thisLoc.SetKmer( kmer );

  return ( otherLoc.atEnd() || thisLoc.atEnd() );
}


bool IsPerfectMatch( KmerPathLoc thisLoc, KmerPathLoc otherLoc ) {
  KmerPathLoc thisCopy(thisLoc), otherCopy(otherLoc);
  return ( thisLoc.GetSegment().Overlaps( otherLoc.GetSegment() ) &&
	   IsPerfectMatchToLeft( thisLoc, otherLoc ) &&
	   IsPerfectMatchToRight( thisCopy, otherCopy ) );
}


bool HavePerfectMatch( const KmerPath& path1, const KmerPath& path2 ) {

  longlong kmer1 = path1.Start(0), kmer2 = path2.Start(0);

  // Check whether path1's first kmer matches in path2:
  for( int seg=0; seg < path2.NSegments(); seg++ ) {
    KmerPathLoc begin1(path1, 0);
    KmerPathLoc loc2(path2, seg);
    if( path2.Segment(seg).Contains(kmer1) && 
        IsPerfectMatchToRight( begin1, loc2 ) )
      return true;
  }
  // Check whether path2's first kmer matches in path1:
  for( int seg=0; seg < path1.NSegments(); seg++ ) {
    KmerPathLoc begin2(path2, 0);
    KmerPathLoc loc1(path1, seg);
    if( path1.Segment(seg).Contains(kmer2) &&
	IsPerfectMatchToRight( begin2, loc1 ) )
      return true;
  }
  return false;
}


Bool ProperOverlapExt( const KmerPath& p1, const KmerPath& p2, int ind1, int ind2 )
{    for ( int i = 0; ; i++ )
     {    if ( i > 0 && p1.Start(ind1+i) != p2.Start(ind2+i) ) return False;
          Bool end = False;
          if ( ind1+i == p1.NSegments( ) - 1 )
          {    if ( p1.Stop(ind1+i) <= p2.Stop(ind2+i) ) break;
               else if ( ind2+i < p2.NSegments( ) - 1 ) return False;
               end = True;    }
          if ( ind2+i == p2.NSegments( ) - 1 )
          {    if ( p2.Stop(ind2+i) <= p1.Stop(ind1+i) ) break;
               else if ( ind1+i < p1.NSegments( ) - 1 ) return False;
               end = True;    }
          if ( !end && p1.Stop(ind1+i) != p2.Stop(ind2+i) ) return False;    }
     for ( int i = 0; ; i++ )
     {    if ( i > 0 && p1.Stop(ind1-i) != p2.Stop(ind2-i) ) return False;
          Bool end = False;
          if ( ind1-i == 0 ) 
          {    if ( p1.Start(ind1-i) >= p2.Start(ind2-i) ) break;
               else if ( ind2-i > 0 ) return False;
               end = True;    }
          if ( ind2-i == 0 )
          {    if ( p2.Start(ind2-i) >= p1.Start(ind1-i) ) break;
               else if ( ind1-i > 0 ) return False;
               end = True;    }
          if ( !end && p1.Start(ind1-i) != p2.Start(ind2-i) ) return False;    }
     return True;    }

template<class TAG> void CreateDatabase(
     const vecKmerPath& paths, vec<TAG>& pathsdb )
{    pathsdb.clear( );
     pathsdb.reserve( paths.sumSizes() );
     for ( size_t i = 0; i < paths.size( ); i++ )
          paths[i].AppendToDatabase( pathsdb, i );
     Prepare(pathsdb);    }

template<class TAG> void CreateDatabase( 
     const vecKmerPath& paths, const vecKmerPath& paths_rc, vec<TAG>& pathsdb )
{    pathsdb.clear( );
     pathsdb.reserve( paths.sumSizes() + paths_rc.sumSizes() );
     for ( size_t i = 0; i < paths.size( ); i++ )
          paths[i].AppendToDatabase( pathsdb, i );
     for ( size_t i = 0; i < paths_rc.size( ); i++ )
          paths_rc[i].AppendToDatabase( pathsdb, ~i ); // ~i == -i-1
     Prepare(pathsdb);    }

template void CreateDatabase( const vecKmerPath& paths, vec<tagged_rpint>& pathsdb );
template void CreateDatabase( const vecKmerPath& paths, 
     const vecKmerPath& paths_rc, vec<tagged_rpint>& pathsdb );
template void CreateDatabase( 
     const vecKmerPath& paths, vec<big_tagged_rpint>& pathsdb );
template void CreateDatabase( const vecKmerPath& paths, 
     const vecKmerPath& paths_rc, vec<big_tagged_rpint>& pathsdb );
template void CreateDatabase( 
     const vecKmerPath& paths, vec<new_tagged_rpint>& pathsdb );
template void CreateDatabase( const vecKmerPath& paths, 
     const vecKmerPath& paths_rc, vec<new_tagged_rpint>& pathsdb );



// MarkReadsWithNovelKmers
// Find which reads contain kmers NOT in this KmerPathDatabase.
template <class TAG>
void
MarkReadsWithNovelKmers( const vecKmerPath & reads, const vec<TAG>& reference,
			 vec<Bool>& answer )
{
  size_t n_reads = reads.size();
  answer.resize( n_reads, False );
  vec<longlong> kmer_locs;
  
  // Loop over all KmerPathIntervals in all reads.
  for ( size_t i = 0; i < n_reads; i++ )
    for ( int j = 0; j < reads[i].NSegments(); j++ ) {
      
      // See whether this KmerPathInterval appears completely in the database.
      // If not, mark this read.
      Contains( reference, reads[i].Segment(j), kmer_locs );
      if ( kmer_locs.empty() ) {
	answer[i] = True;
	break;
      }
    }
}

// Template instantiations.
template void MarkReadsWithNovelKmers( const vecKmerPath & reads, const vec<tagged_rpint>& segs, vec<Bool>& answer );


// MarkReadsWithNovelKmers
// Find which reads contain kmers NOT in this KmerPathDatabase or not present in correctly adjacent unipaths
template <class TAG>
void
MarkReadsWithNovelKmers2( const vecKmerPath & reads, const vec<TAG>& reference, const digraph& adj_graph,
			 vec<Bool>& answer )
{
  size_t n_reads = reads.size();
  answer.resize( n_reads, False );
  vec<longlong> kmer_locs;
  
  // Loop over all KmerPathIntervals in all reads.
  for ( size_t ir = 0; ir < n_reads; ir++ ){
    vec<int> uorder;
    for ( int j = 0; j < reads[ir].NSegments(); j++ ) {
      
      // See whether this KmerPathInterval appears completely in the database.
      // If not, mark this read.
      Contains( reference, reads[ir].Segment(j), kmer_locs );
      if ( kmer_locs.empty() ) {
	answer[ir] = True;
	break;
      }else{
	for ( size_t l = 0; l < kmer_locs.size(); l++ ){
	  const tagged_rpint& t = reference[kmer_locs[l]];
	  int tid = t.PathId();
	  if ( tid < 0 )
	    FatalErr("Unanticipated negative result.");
	  if ( uorder.size() == 0 || uorder.back() != tid )
	    uorder.push_back(tid);
	}
      }
    }
    if ( uorder.size() > 1 ){
      // check this unipath alignment order in the unipath adjacency graph
      for ( size_t k = 0; k < uorder.size() -1; k++ ){
	int uid1 = uorder[k], uid2 = uorder[k+1];
	if ( ! adj_graph.HasEdge( uid1, uid2 ) ){
	  answer[ir] = True;
	  //cout << "\nHERE IT IS: found unipath kmers not adjacent in the unipath adjacency graph: " << endl;
	  break;
	}
      }
    }
    //cout << "uorder: "; uorder.Print(cout); cout << endl;
    //PRINT( ToStringBool(answer[ir]) );
  }
}

// Template instantiations.
template void MarkReadsWithNovelKmers2( const vecKmerPath & reads, const vec<tagged_rpint>& segs, const digraph&, vec<Bool>& answer );


Bool SubContig( const KmerPath& p, const vecKmerPath& paths, 
     const vecKmerPath& paths_rc, const vec<tagged_rpint>& xpathsdb,
     const HashSimple* pPathsToIgnore, read_id_t pathToIgnore )
{    longlong seg = 0, pos_on_seg = 0, last_seg = p.NSegments( ) - 1;
     Bool at_end = False;
     while(1)
     {    kmer_id_t x = p.Segment(seg).Start( ) + pos_on_seg;
          vec<longlong> instances;
	  Contains( xpathsdb, x, instances );
	  instances.ReverseMe(); // this creates the effect of a KmerOccurrenceIterBwd
          Bool succeed = False;
	  for ( size_t i = 0; i < instances.size(); i++ ) {
	       const tagged_rpint & t = xpathsdb[ instances[i] ];
               longlong id2 = t.PathId( );
               read_id_t readId2 = t.ReadId();
               if ( pathToIgnore == readId2  ||
		    pPathsToIgnore &&
                    pPathsToIgnore->Has( readId2 ) )
                 continue;

               const KmerPath& q = ( id2 >= 0 ? paths[id2] : paths_rc[-id2-1] );
               longlong qseg = t.PathPos( );
               if ( !ProperOverlapExt( p, q, seg, qseg ) ) continue;
               int extension = q.Segment(qseg).Stop( ) - x;
               for ( int u = qseg + 1; u < q.NSegments( ); u++ )
                    extension += q.Segment(u).Length( );
               if ( extension > 0 ) succeed = True;
               while( extension > 0 )
               {    int last_pos = p.Segment(seg).Length( ) - 1;
                    if ( pos_on_seg < last_pos )
                    {    if ( last_pos - pos_on_seg <= extension )
                         {    extension -= ( last_pos - pos_on_seg );
                              pos_on_seg = last_pos;    }
                         else
                         {    pos_on_seg += extension;
                              extension = 0;    }    }
                    else
                    {    if ( seg == last_seg ) 
                         {    at_end = True;
                              break;    }
                         ++seg;
                         pos_on_seg = 0;
                         --extension;    }    }
               if (succeed) break;    }
          if ( !succeed || at_end ) break;    }
     return at_end;    }

longlong KmerPath::GetKmer( int n ) const
{    int seg = 0;
     while(1)
     {    if ( n < Length(seg) ) return Start(seg) + n;
          n -= Length(seg);
          ++seg;
          ForceAssertLt( seg, NSegments( ) );    }    }


vecKmerPathIndex::vecKmerPathIndex( const vecKmerPath& _gpaths, const vecKmerPath& _gpaths_rc ):
  gpaths(_gpaths), gpaths_rc(_gpaths_rc), segmentStarts(_gpaths.size()), segmentStarts_rc(_gpaths.size()) {
  
  ForceAssert( _gpaths.size() == _gpaths_rc.size() );
  
  // For each path, build an index of where the segments of this path start.
  for ( int rc = 0; rc < 2; rc++ ) {
    const vecKmerPath& paths = rc ? gpaths : gpaths_rc;
    vec< vec< genome_part_pos_t > >& segStarts = rc ? segmentStarts : segmentStarts_rc;
    for ( size_t i = 0; i < paths.size( ); i++ )
      {    segStarts[i].reserve( paths[i].NSegments( ) + 1 );
      segStarts[i].push_back(0);
      for ( int j = 0; j < paths[i].NSegments( ); j++ )
	{
	  ForceAssert( paths[i].Segment(j).Length( ) > 0 );
	  segStarts[i].push_back( segStarts[i].back( ) 
				     + paths[i].Segment(j).Length( ) );    }    }
  }
    
}

// Method: FindLoc
// 
// Given a position in the path, find the segment and the location within the segment
// for the kmer starting at that position.
KmerPathLoc vecKmerPathIndex::FindLoc( longlong pathId, int posInPath, Bool orient ) const {
  const KmerPath& path = (orient == ORIENT_FW ? gpaths : gpaths_rc)[pathId];
  const vec< genome_part_pos_t > & segStarts =
    (orient == ORIENT_FW  ? segmentStarts : segmentStarts_rc)[pathId];

  cout << " looking up loc " << posInPath << endl;
  
  // segStarts[j] contains the start of segment j.
  // find the first segment in path, such that posInPath falls on its beginning or to the right of its
  // beginning.
  vec< genome_part_pos_t >::const_iterator
    segContainingPos = upper_bound( segStarts.begin(),
				    segStarts.end(),
				    posInPath );

  int segId = segContainingPos - segStarts.begin() - 1;
  int segOffset = posInPath - segStarts[segId];
  ForceAssert( &path == &((orient == ORIENT_FW ? gpaths : gpaths_rc)[pathId]) );
  return KmerPathLoc( path, segId, segOffset);
}

template void digraphE<KmerPath>::GiveEdgeNewFromVx(int, int, int);
template void digraphE<KmerPath>::GiveEdgeNewToVx(int, int, int);
template Bool digraphE<KmerPath>::EdgePaths(int, int, vec<vec<int> >&, const int,
     const int, const int);
template vec<int> digraphE<KmerPath>::EdgesBetween(int, int) const;
template int digraphE<KmerPath>::EdgeObjectIndexToToIndex(int, int) const;
template void digraphE<KmerPath>::DeleteEdges(vec<int> const&);
template int digraphE<KmerPath>::EdgeObjectIndexToFromIndex(int, int) const;
template int digraphE<KmerPath>::EdgeObjectIndexByIndexTo(int, int) const;
template const KmerPath& digraphE<KmerPath>::EdgeObject(int) const;
template KmerPath& digraphE<KmerPath>::EdgeObjectMutable(int);
template void digraphE<KmerPath>::Clear( );
template void digraphE<KmerPath>::AddVertices(int);

template class SmallVec< KmerPathInterval, MempoolAllocator<KmerPathInterval> >;
template class OuterVec<KmerPath>;

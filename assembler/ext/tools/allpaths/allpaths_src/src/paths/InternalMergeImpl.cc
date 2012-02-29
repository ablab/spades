///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <set>

#include "Basevector.h"
#include "Equiv.h"
#include "Intvector.h"
#include "TaskTimer.h"
#include "math/Functions.h"
#include "paths/ImproperMerge.h"
#include "paths/InternalMergeImpl.h"
#include "graph/DigraphTemplate.h"




/*******************************************************************************
 *
 * Classes used ONLY in this module:
 *
 * KmerPathLocAlt
 * MultiKmerPathLoc
 * PointedSubPathPair
 * KmerPathSplitPoint
 * KmerPathCorrespondence
 *
 ******************************************************************************/


// A KmerPathLocAlt stores the position of a kmer on a KmerPath as three numbers:
// the index of the KmerPathInterval in the KmerPath, the index of the kmer
// in the KmerPathInterval, and the index of the kmer in the entire KmerPath
// (ignoring gaps).  The third item is redundant if you know the KmerPath.  This
// class should not be used for storing positions in gaps.
//
// The second constructor takes a KmerPath as one of its inputs and uses it to
// calculate the index of the kmer in the KmerPath.  The third constructor is
// from a KmerPathLoc object.

class KmerPathLocAlt {

     public:

     KmerPathLocAlt( ) : interval_(0), pos_on_interval_(0), pos_on_path_(0) { }
     KmerPathLocAlt( const KmerPath& p, int interval, int pos_on_interval )
          : interval_(interval), pos_on_interval_(pos_on_interval)
     {    pos_on_path_ = pos_on_interval;
          for ( int i = 0; i < interval; i++ )
          {    if ( p.isGap(i) ) continue;
               pos_on_path_ += p.Length(i);    }    }
     KmerPathLocAlt( const KmerPath& p, int interval, int pos_on_interval,
          int pos_on_path ) : interval_(interval),
          pos_on_interval_(pos_on_interval), pos_on_path_(pos_on_path) { }
     KmerPathLocAlt( const KmerPathLoc& l )
          : interval_( l.GetIndex( ) ),
          pos_on_interval_( l.GetLoc( ) )
     {    ForceAssert( !l.isGap( ) );
          const KmerPath& p = l.GetPath( );
          pos_on_path_ = pos_on_interval_;
          for ( int i = 0; i < interval_; i++ )
          {    if ( p.isGap(i) ) continue;
               pos_on_path_ += p.Length(i);    }    }

     int Interval( ) const { return interval_; }
     int PosOnInterval( ) const { return pos_on_interval_; }
     int PosOnPath( ) const { return pos_on_path_; }

     void Reverse( const KmerPath& p ) {
       pos_on_interval_ = p.Segment(interval_).Length( ) - pos_on_interval_ - 1;
       interval_ = p.NSegments( ) - interval_ - 1;
       pos_on_path_ = pos_on_interval_;
       for ( int i = 0; i < interval_; i++ ) {
	 if ( p.isGap(i) ) continue;
	 pos_on_path_ += p.Length(i);
       }
     }

     friend Bool operator<( const KmerPathLocAlt& l1, const KmerPathLocAlt& l2 )
     {    if ( l1.Interval( ) < l2.Interval( ) ) return True;
          if ( l1.Interval( ) > l2.Interval( ) ) return False;
          return l1.PosOnInterval( ) < l2.PosOnInterval( );    }
     friend Bool operator>( const KmerPathLocAlt& l1, const KmerPathLocAlt& l2 )
     {    return l2 < l1;    }

     private:

     int interval_, pos_on_interval_, pos_on_path_;

};

TRIVIALLY_SERIALIZABLE(KmerPathLocAlt);







// A MultiKmerPathLoc represents a position on a vecKmerPath.
class MultiKmerPathLoc {

public:

  // Boring constructors.

  MultiKmerPathLoc( ) { }
  MultiKmerPathLoc( int edge, int segment_on_edge, int kmer_on_segment )
    : edge_(edge), segment_on_edge_(segment_on_edge),
      kmer_on_segment_(kmer_on_segment) { }
  MultiKmerPathLoc( int edge, const KmerPathLocAlt& l )
    : edge_(edge), segment_on_edge_( l.Interval( ) ),
      kmer_on_segment_( l.PosOnInterval( ) ) { }

  // The following constructor takes as input a vecKmerPath and a position
  // on its merger (which is a KmerPath), and from that generates a position
  // on the vecKmerPath.

  MultiKmerPathLoc( const vecKmerPath& m, const KmerPathLoc& l );

  KmerPathLoc PosOnMerger( const vecKmerPath& m, const KmerPath& merger ) const;

  longlong Kmer( const vecKmerPath& m )
  {    return m[ Edge( ) ].Segment( SegmentOnEdge( ) ).Start( ) 
      + KmerOnSegment( );    }

  // Access the data:

  int Edge( ) const { return edge_; }
  int SegmentOnEdge( ) const { return segment_on_edge_; }
  int KmerOnSegment( ) const { return kmer_on_segment_; }

  // Set data.

  void SetEdge( int e ) { edge_ = e; }
  void SetSegmentOnEdge( int s ) { segment_on_edge_ = s; }
  void SetKmerOnSegment( int k ) { kmer_on_segment_ = k; }

  void AddToSegmentOnEdge( int a ) { segment_on_edge_ += a; }

  // Test for ends.  The right end requires more data because this class doesn't
  // own enough information.

  Bool LeftEnd( ) const
  {    return Edge( ) == 0 && SegmentOnEdge( ) == 0 && KmerOnSegment( ) == 0;    }

  Bool RightEnd( const EmbeddedSubPath<KmerPath>& p ) const
  {    const KmerPath& kp = p.EdgeObject( Edge( ) );
    int lastseg = kp.NSegments( ) - 1;
    return Edge( ) == p.NEdges( ) - 1 && SegmentOnEdge( ) == lastseg
      && KmerOnSegment( ) == kp.Segment(lastseg).Length( ) - 1;    }

  Bool AtRightEndOfEdge( const vecKmerPath& m ) const
  {    const KmerPath& km = m[ Edge( ) ];
    int lastseg = km.NSegments( ) - 1;
    return SegmentOnEdge( ) == lastseg
      && KmerOnSegment( ) == km.Segment(lastseg).Length( ) - 1;    }

  friend Bool operator<=(const MultiKmerPathLoc& l1, const MultiKmerPathLoc& l2);
  friend Bool operator>=( const MultiKmerPathLoc& l1, const MultiKmerPathLoc& l2 )
  {    return l2 <= l1;    }

  friend Bool operator==( const MultiKmerPathLoc& l1, const MultiKmerPathLoc& l2 )
  {    return l1.Edge( ) == l2.Edge( ) 
      && l1.SegmentOnEdge( ) == l2.SegmentOnEdge( )
      && l1.KmerOnSegment( ) == l2.KmerOnSegment( );    }

  friend ostream& operator<<( ostream& out, const MultiKmerPathLoc& l )
  {    return out << l.Edge( ) << "." << l.SegmentOnEdge( ) << "."
		  << l.KmerOnSegment( );    }

private:

  int edge_;
  int segment_on_edge_;
  int kmer_on_segment_;

};


MultiKmerPathLoc::MultiKmerPathLoc( const vecKmerPath& m, const KmerPathLoc& l )
{    edge_ = 0;
     segment_on_edge_ = l.GetIndex( );
     kmer_on_segment_ = l.GetLoc( );
     for ( size_t u = 0; u < m.size( ); u++ )
     {    int segs_in_edge = m[u].NSegments( );
          if ( segs_in_edge > segment_on_edge_ ) break;
          segment_on_edge_ -= segs_in_edge;
          ++edge_;    }    }

KmerPathLoc MultiKmerPathLoc::PosOnMerger(
     const vecKmerPath& m, const KmerPath& merger ) const
{    int segment = SegmentOnEdge( );
     for ( int i = 0; i < Edge( ); i++ )
          segment += m[i].NSegments( );
     return KmerPathLoc( merger, segment, KmerOnSegment( ) );    }


Bool operator<=( const MultiKmerPathLoc& l1, const MultiKmerPathLoc& l2 )
{    if ( l1.Edge( ) < l2.Edge( ) ) return True;
     if ( l1.Edge( ) > l2.Edge( ) ) return False;
     if ( l1.SegmentOnEdge( ) < l2.SegmentOnEdge( ) ) return True;
     if ( l1.SegmentOnEdge( ) > l2.SegmentOnEdge( ) ) return False;
     if ( l1.KmerOnSegment( ) < l2.KmerOnSegment( ) ) return True;
     if ( l1.KmerOnSegment( ) > l2.KmerOnSegment( ) ) return False;
     return True;    }




// A PointedSubPathPair consists of two EmbeddedSubPaths of a HyperKmerPath,
// and positions on them, referring to equal kmers.

class PointedSubPathPair {

public:

  void TestValid( ) const; // incomplete test

  PointedSubPathPair( ) { }
  PointedSubPathPair( const EmbeddedSubPath<KmerPath>& a1, 
		      const EmbeddedSubPath<KmerPath>& a2, const MultiKmerPathLoc& l1, 
		      const MultiKmerPathLoc& l2 )
    : a1_(a1), a2_(a2), l1_(l1), l2_(l2)
  {    TestValid( );    }

  const EmbeddedSubPath<KmerPath>& Path1( ) const { return a1_; }
  const EmbeddedSubPath<KmerPath>& Path2( ) const { return a2_; }
  const EmbeddedSubPath<KmerPath>& Path( int i ) const 
  {    return ( i == 0 ? a1_ : a2_ );    }
  EmbeddedSubPath<KmerPath>& PathMutable( int i )
  {    return ( i == 0 ? a1_ : a2_ );    }

  const MultiKmerPathLoc& Loc1( ) const { return l1_; }
  const MultiKmerPathLoc& Loc2( ) const { return l2_; }
  const MultiKmerPathLoc& Loc( int i ) const
  {    return ( i == 0 ? l1_ : l2_ );    }
  MultiKmerPathLoc& LocMutable( int i )
  {    return ( i == 0 ? l1_ : l2_ );    }

  void PrintSummary( ostream& out, int i, int K ) const;
  void PrintSummary( ostream& out, int K ) const
  {    PrintSummary( out, 0, K );
    cout << "\n";
    PrintSummary( out, 1, K );    
    cout << "\n";    }

  friend bool operator== ( const PointedSubPathPair& lhs, const PointedSubPathPair& rhs ) {
    return ( lhs.a1_ == rhs.a1_ &&
	     lhs.a2_ == rhs.a2_ &&
	     lhs.l1_ == rhs.l1_ &&
	     lhs.l2_ == rhs.l2_ );
  }

private:

  EmbeddedSubPath<KmerPath> a1_, a2_;
  MultiKmerPathLoc l1_, l2_;

};


void PointedSubPathPair::TestValid( ) const
{    const KmerPath& p1 = Path1( ).EdgeObject( Loc1( ).Edge( ) );
     const KmerPath& p2 = Path2( ).EdgeObject( Loc2( ).Edge( ) );
     int seg1 = Loc1( ).SegmentOnEdge( ), seg2 = Loc2( ).SegmentOnEdge( );
     ForceAssertEq(
          p1.Segment(seg1).Start( ) + Loc1( ).KmerOnSegment( ),
          p2.Segment(seg2).Start( ) + Loc2( ).KmerOnSegment( ) );    }

void PointedSubPathPair::PrintSummary( ostream& out, int i, int K ) const
{    const EmbeddedSubPath<KmerPath>& e = Path(i);
     for ( int j = 0; j < e.NVertices( ); j++ )
     {    out << e.Vertex(j);
          if ( j == 0 && e.IsFromSource() ) out << "[source]";
          if ( j < e.NVertices( ) - 1 )
          {    const KmerPath& p = e.EdgeObject(j);
               float len = 0, dev = 0;
               for ( int u = 0; u < p.NSegments( ); u++ )
               {    const KmerPathInterval& x = p.Segment(u);
                    if ( x.isSeq( ) ) len += x.Length( );
                    else
                    {    len += float( x.Maximum( ) + x.Minimum( ) ) / 2.0;
                         dev += float( x.Maximum( ) - x.Minimum( ) ) / 2.0;    }    }
               out << " ---(" << len;
               if ( dev > 0 ) out << "~" << dev;
               out << ")--> ";    }
          else if ( e.IsToSink() ) out << "[sink]";
     }
}





// A KmerPathSplitPoint represents a point in a KmerPath at which it could be split
// into two pieces (allowing for the possibility that one is null).  We do not allow
// splitting within a gap.
//
// The data is tracked via counting the number of complete kmer intervals to the 
// left of the split point, plus zero or more kmers in the partial kmer interval
// to the immediate left of the split point.

class KmerPathSplitPoint {

public:

  KmerPathSplitPoint( ) 
    : intervals_to_left_(-1), kmers_in_partial_to_left_(-1) { }
  KmerPathSplitPoint( int intervals_to_left, int kmers_in_partial_to_left )
    : intervals_to_left_(intervals_to_left), 
      kmers_in_partial_to_left_(kmers_in_partial_to_left) { }

  Bool Initialized( ) const { return intervals_to_left_ >= 0; }

  int IntervalsToLeft( ) const { return intervals_to_left_; }
  int KmersInPartialToLeft( ) const { return kmers_in_partial_to_left_; }

  void SetIntervalsToLeft( int x ) { intervals_to_left_ = x; }

  friend Bool operator<( const KmerPathSplitPoint& s1, 
			 const KmerPathSplitPoint& s2 )
  {    if ( s1.IntervalsToLeft( ) < s2.IntervalsToLeft( ) ) return True;
    if ( s1.IntervalsToLeft( ) > s2.IntervalsToLeft( ) ) return False;
    if ( s1.KmersInPartialToLeft( ) < s2.KmersInPartialToLeft( ) ) return True;
    return False;    }

  friend Bool operator==( const KmerPathSplitPoint& s1, 
			  const KmerPathSplitPoint& s2 )
  {    return s1.IntervalsToLeft( ) == s2.IntervalsToLeft( )
      && s1.KmersInPartialToLeft( ) == s2.KmersInPartialToLeft( );    }

  // Between: return the part of a KmerPath between two split points.  We don't
  // allow the split points to be equal.

  friend KmerPath Between( const KmerPath& p, const KmerPathSplitPoint& left,
			   const KmerPathSplitPoint& right );

  friend ostream& operator<<( ostream& out, const KmerPathSplitPoint& s )
  {    return out << s.IntervalsToLeft( ) << "." 
		  << s.KmersInPartialToLeft( );    }

private:

  int intervals_to_left_;
  int kmers_in_partial_to_left_;

};


// Class: KmerPathCorrespondence
//
// A KmerPathCorrespondence is just the index of one KmerPath, and the position
// of a kmer on it, and the index of another KmerPath, and the position of a kmer
// on it, such that the kmers are the same.
//
// In the constructor, the edge on the first path goes from vertex from1 to vertex
// to1, and the edge on the second path goes from vertex from2 to vertex to2.  The
// indices of the edges in digraphE::edges_ are id1 and id2.  The index of the first
// edge in from_[from1] is from_index1 and the index of the second edge in
// from_[from2] is from_index2.
//
// Typical usage: a KmerPathCorrespondence between two edges of a <HyperKmerPath>
// indicates that the two edges may potentially be merged to create a more compact
// representation of the closure set represented by the HyperKmerPath.
class KmerPathCorrespondence {

     public:

     KmerPathCorrespondence( ) { }
     KmerPathCorrespondence( int id1, int from1, int from_index1, int to1, 
          const KmerPathLocAlt& l1, int id2, int from2, int from_index2, int to2, 
          const KmerPathLocAlt& l2 ) 
          : id1_(id1), id2_(id2), from1_(from1), from2_(from2), 
          from_index1_(from_index1), from_index2_(from_index2), to1_(to1), 
          to2_(to2), l1_(l1), l2_(l2) { }

     void Swap( );

     int Id1( ) const { return id1_; }
     int Id2( ) const { return id2_; }
     int From1( ) const { return from1_; }
     int From2( ) const { return from2_; }
     int FromIndex1( ) const { return from_index1_; }
     int FromIndex2( ) const { return from_index2_; }
     int To1( ) const { return to1_; }
     int To2( ) const { return to2_; }
     KmerPathLocAlt Pos1( ) const { return l1_; }
     KmerPathLocAlt Pos2( ) const { return l2_; }
     int Offset( ) const { return l1_.PosOnPath( ) - l2_.PosOnPath( ); }

     private:
     
     int id1_, id2_;
     int from1_, from2_;
     int from_index1_, from_index2_;
     int to1_, to2_;
     KmerPathLocAlt l1_, l2_;

};

void KmerPathCorrespondence::Swap( )
{    swap( id1_, id2_ );
     swap( from1_, from2_ );
     swap( from_index1_, from_index2_ );
     swap( to1_, to2_ );
     swap( l1_, l2_ );    }










// SharedKmers: find all pairs of edges i <= j which share a kmer, returning
// a vec<KmerPathCorrespondence> to capture this info.  We do not
// return multiple position pairs which are simple shifts of each other along
// path segments.  Also, we do not return identity matches.
void SharedKmers( const HyperKmerPath & hkp, vec<KmerPathCorrespondence>& shares,
		  int min_improper )
{    shares.clear( );

     // Build path database.

     vec<Bool> used;
     int n_edges = hkp.EdgeObjectCount( );
     hkp.Used(used);
     int nused = Sum(used);
     vec<int> from( n_edges ), from_index( n_edges ), to( n_edges );
     for ( int v = 0; v < hkp.N( ); v++ )
     {    for ( int j = 0; j < hkp.From(v).isize( ); j++ )
          {    int w = hkp.From(v)[j];
               int e = hkp.EdgeObjectIndexByIndexFrom( v, j );
               from[e] = v;
               from_index[e] = j;
               to[e] = w;    }    }
     vec<tagged_rpint> segs;
     segs.reserve(nused);
     for ( int i = 0; i < n_edges; i++ )
          if ( used[i] ) hkp.EdgeObject(i).AppendToDatabase( segs, i );
     Prepare(segs);


     vec<int> edgelength( n_edges );
     int n_kpis = 0;
     for ( int i = 0; i < n_edges; i++ ) {
       edgelength[i] = hkp.EdgeLength( i );
       n_kpis += hkp.EdgeObject( i ).NSegments( );
     }
     
     // Make a table of the starting locations of each KmerPathInterval on each KmerPath.
     VecIntVec starts;
     starts.reserve( n_edges );
     IntVec kp_starts;
     for ( int i = 0; i < n_edges; i++ ) {
       kp_starts.clear( );
       int start = 0;
       for ( int j = 0; j < hkp.EdgeObject( i ).NSegments( ); j++ ) {
	 kp_starts.push_back( start );
	 start += hkp.EdgeObject( i ).Length( j );
       }
       starts.push_back( kp_starts );
     }
     vec<int> share_lengths;
     
     // Search path database for shared kmers.
     for ( int i = 0; i < segs.isize( ); i++ )
     {    int id1 = segs[i].PathId( ), pos1 = segs[i].PathPos( );
          vec<longlong> matches;
          const KmerPathInterval& rpi1 = hkp.EdgeObject(id1).Segment(pos1);
          Contains( segs, rpi1, matches );
          for ( int j = 0; j < matches.isize( ); j++ )
          {    const tagged_rpint& t = segs[ matches[j] ];
               int id2 = t.PathId( ), ppos2 = t.PathPos( );
               if ( id2 < id1 ) continue;
               if ( id1 == id2 && pos1 == ppos2 ) continue;
               const KmerPathInterval& rpi2 = hkp.EdgeObject(id2).Segment(ppos2);

               // Now we have two overlapping path intervals.  Exhibit positions
               // on each which correspond to the same kmer.

               int ind1, ind2;
               longlong start1 = rpi1.Start( ), stop1 = rpi1.Stop( );
               longlong start2 = rpi2.Start( ), stop2 = rpi2.Stop( );
               if ( start1 >= start2 )
               {    ind1 = 0, ind2 = int(start1 - start2);    }
               else
               {    ind2 = 0, ind1 = int(start2 - start1);    }

               // Reject if the kmers to the left match.

               if ( ind1 == 0 && ind2 == 0 && pos1 > 0 && ppos2 > 0 )
               {    if ( hkp.EdgeObject(id1).Stop(pos1-1)
                         == hkp.EdgeObject(id2).Stop(ppos2-1) )
                    {    continue;    }    }

               // Generate the share.

               int pathpos1 = starts[id1][pos1] + ind1;
               int pathpos2 = starts[id2][ppos2] + ind2;
               KmerPathLocAlt l1( hkp.EdgeObject(id1), pos1, ind1, pathpos1 );
               KmerPathLocAlt l2( hkp.EdgeObject(id2), ppos2, ind2, pathpos2 );
               KmerPathCorrespondence s( id1, from[id1], from_index[id1],
                    to[id1], l1, id2, from[id2], from_index[id2], to[id2], l2 );

               // Remove shares whose perfect extensions P are improper and such
               // that length(P) is less than min_improper.

               int e1 = s.Id1( ), e2 = s.Id2( );
               KmerPathLoc scanStart1( hkp.EdgeObject(e1), s.Pos1( ).Interval( ),
                    s.Pos1( ).PosOnInterval( ) );
               KmerPathLoc scanStart2( hkp.EdgeObject(e2), s.Pos2( ).Interval( ),
                    s.Pos2( ).PosOnInterval( ) );

               KmerPathLoc rightScan1( scanStart1 ), rightScan2( scanStart2 );
               ScanRightPerfectMatch( rightScan1, rightScan2 );
               int right_ext = KmersInInterval( scanStart1, rightScan1 ) - 1;

               KmerPathLoc leftScan1( scanStart1 ), leftScan2( scanStart2 );
               ScanLeftPerfectMatch( leftScan1, leftScan2 );
               int left_ext = KmersInInterval( leftScan1, scanStart1 ) - 1;

               int pos1 = s.Pos1( ).PosOnPath( ), pos2 = s.Pos2( ).PosOnPath( );
               int len1 = edgelength[e1], len2 = edgelength[e2];
               int left1 = pos1 - left_ext, right1 = pos1 + right_ext + 1;
               int left2 = pos2 - left_ext, right2 = pos2 + right_ext + 1;
               Bool left_proper = ( left1 == 0 || left2 == 0 );
               Bool right_proper = ( right1 == len1 || right2 == len2 );
               Bool proper = ( left_proper && right_proper );
	       int share_length = left_ext + right_ext + 1;
               if ( !proper && share_length < min_improper ) continue;

               // Save the share.

               shares.push_back(s);
	       share_lengths.push_back( share_length );
	  }
     }

     // Sort shares.

     SortSync( share_lengths, shares );
}



// BAD BAD BAD:

int MinLen( const KmerPathInterval& I )
{    if ( I.isSeq( ) ) return I.Length( );
     else return 1000000;    }  // VERY SLOPPY.

// Advance: move a KmerPathLoc to the next kmer.  The starting location must not be
// in a gap.

void Advance( KmerPathLoc& l )
{    const KmerPath& p = l.GetPath( );
     ForceAssert( !p.Segment( l.GetIndex( ) ).isGap( ) );
     if ( l.GetLoc( ) + 1 < p.Segment( l.GetIndex( ) ).Length( ) )
          l.SetLoc( l.GetLoc( ) + 1 );
     else
     {    ForceAssertLt( l.GetIndex( ) + 1, p.NSegments( ) );
          l.SetIndex( l.GetIndex( ) + 1 );
          l.SetLoc(0);    }    }

Bool InRange( const KmerPathLoc& l )
{    if ( l.GetIndex( ) < 0 ) return False;
     if ( l.GetIndex( ) >= l.GetPath( ).NSegments( ) ) return False;
     if ( l.GetLoc( ) < 0 ) return False;
     if ( l.GetLoc( ) >= l.GetPath( ).Segment( l.GetIndex( ) ).Length( ) )
          return False;
     return True;    }


// Create a vecKmerPath out of an EmbeddedSubPath of KmerPaths.
void
CondenseToVecKmerPath( const EmbeddedSubPath<KmerPath>& in, vecKmerPath & out )
{
  // Pre-reserve memory.
  int n_edges = in.NEdges( ), n_kpis = 0;
  for ( int i = 0; i < n_edges; i++ )
    n_kpis += in.EdgeObject( i ).NSegments( );
  out.Reserve( n_edges + 2*n_kpis, n_edges ); // each KmerPathInterval takes up 8 bytes, or 2 ints
  
  // Create the vecKmerPath.
  for ( int i = 0; i < n_edges; i++ )
    out.push_back( in.EdgeObject(i) );
}

// Concatenate the KmerPaths in a vecKmerPath.  We deliberately avoid
// concatenating abutting KmerPathIntervals because that would make it
// harder to do accounting when using vecKmerPaths.
KmerPath MergeKmerPaths( const vecKmerPath& in )
{
  KmerPath out;
  // Reserve memory for the merged KmerPath.
  int size = 0;
  for ( size_t i = 0; i < in.size( ); i++ )
    size += in[i].NSegments( );
  out.Reserve( size );
  // Create the merged KmerPath.
  for ( size_t i = 0; i < in.size( ); i++ )
    for ( int j = 0; j < in[i].NSegments( ); j++ )
      out.AddSegmentNoConcatenate( in[i].Segment(j) );
  return out;
}




// AlignSubpaths.  Given a PointedSubPathPair, look for partial alignments 
// extending the matching kmers, and return the number naligns of alignments 
// found.  If naligns = 1, return details as follows:
// - ends[0] and ends[1]: positions of leftmost and rightmost kmers of 
//   partial alignment on first and second sequences;
// - M: merged path
// - i0,...,im: positions of first sequence's vertices on M, if in
//   the aligning part;
// - j0,...,jn: positions of second sequence's vertices on M, if in
//   the aligning part.
// CHANGED: if n > 1, return details for the FIRST alignment found.

// Note that if AlignSubpaths finds more than one path, it returns the first one.
// this can't be optimal.
int AlignSubpaths( const PointedSubPathPair& p,
     vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> >& ends,
     KmerPath& M, vec<KmerPathSplitPoint>& i, vec<KmerPathSplitPoint>& j,
     const NegativeGapValidator* ngv, bool calculate_ij = false )
{
     i.clear( ), j.clear( );

     // Check inputs. (commented out to save runtime)

     //ForceAssert( p.Path1( ).NEdges( ) > 0 && p.Path2( ).NEdges( ) > 0 );
     //ForceAssert( p.Path1( ).EdgeObject( 0 ).NSegments( ) > 0 );
     //ForceAssert( p.Path2( ).EdgeObject( 0 ).NSegments( ) > 0 );

     // Test for equality.  (Don't deal with end gap case.)

     if ( p.Path1( ) == p.Path2( ) && p.Loc1( ) == p.Loc2( ) &&
	  !p.Path1( ).EdgeObject( p.Path1( ).NEdges( ) ).LastSegment( ).isGap( ) )
     {
          vecKmerPath m1;
          CondenseToVecKmerPath( p.Path1( ), m1 );
          M = MergeKmerPaths( m1 );
          ends.resize(2);
          ends[0].first = ends[1].first = MultiKmerPathLoc( 0, 0, 0 );
          ends[0].second = ends[1].second
               = MultiKmerPathLoc( m1.size( ) - 1, m1.back( ).NSegments( ) - 1,
                    m1.back( ).LastSegment( ).Length( ) - 1 );
          i.resize( p.Path1( ).NVertices( ) ); j.resize( p.Path2( ).NVertices( ) );
          int sum = 0;
          for ( int u = 0; u < p.Path1( ).NVertices( ); u++ )
          {    i[u] = j[u] = KmerPathSplitPoint( sum, 0 );
               if ( u < p.Path1( ).NVertices( ) - 1 )
                    sum += m1[u].NSegments( );    }
          return 1;
     }

     // Define KmerPaths to be aligned.
     vecKmerPath m1, m2;
     CondenseToVecKmerPath( p.Path1( ), m1 );
     CondenseToVecKmerPath( p.Path2( ), m2 );

     // Create merged paths.
     KmerPath p1 = MergeKmerPaths( m1 ), p2 = MergeKmerPaths( m2 );

     // Check for two gaps in a row.  For now, reject such cases.  I don't know how
     // often this happens or even if it does happen.

     for ( int u = 0; u < p1.NSegments( ) - 1; u++ )
       if ( p1.Segment(u).isGap( ) && p1.Segment(u+1).isGap( ) )
	 return 0;
     for ( int u = 0; u < p2.NSegments( ) - 1; u++ )
       if ( p2.Segment(u).isGap( ) && p2.Segment(u+1).isGap( ) )
	 return 0;

     // Define positions of matching kmers on p1 and p2.
     KmerPathLoc w1 = p.Loc1( ).PosOnMerger( m1, p1 );
     KmerPathLoc w2 = p.Loc2( ).PosOnMerger( m2, p2 );
     ForceAssertEq( w1.GetKmer( ), w2.GetKmer( ) );

     // If there is a gap at the beginning or end of either m1 or m2, remove the
     // gaps now, and put them back at the end of AlignSubpaths.  We do this
     // because ImproperMerge doesn't handle nonstandard paths.
     KmerPathInterval front_gap1(0,0), front_gap2(0,0);
     KmerPathInterval back_gap1(0,0), back_gap2(0,0);
     Bool have_front_gap1 = m1[0].Segment(0).isGap( );
     Bool have_front_gap2 = m2[0].Segment(0).isGap( );
     Bool have_back_gap1 = m1.back( ).LastSegment( ).isGap( );
     Bool have_back_gap2 = m2.back( ).LastSegment( ).isGap( );
     if (have_front_gap1)
     {    front_gap1 = m1[0].Segment(0);
          KmerPath m1tail;
          for ( int u = 1; u < m1[0].NSegments( ); u++ )
               m1tail.AddSegmentNoConcatenate( m1[0].Segment(u) );
          m1[0] = m1tail;
          p1 = MergeKmerPaths( m1 );
          w1 = KmerPathLoc( p1, w1.GetIndex( ) - 1, w1.GetLoc( ) );    }
     if (have_front_gap2)
     {    front_gap2 = m2[0].Segment(0);
          KmerPath m2tail;
          for ( int u = 1; u < m2[0].NSegments( ); u++ )
               m2tail.AddSegmentNoConcatenate( m2[0].Segment(u) );
          m2[0] = m2tail;
          p2 = MergeKmerPaths( m2 );
          w2 = KmerPathLoc( p2, w2.GetIndex( ) - 1, w2.GetLoc( ) );    }
     if (have_back_gap1)
     {    back_gap1 = m1.back( ).LastSegment( );
          KmerPath m1head;
          for ( int u = 0; u < m1.back( ).NSegments( ) - 1; u++ )
               m1head.AddSegmentNoConcatenate( m1.back( ).Segment(u) );
          m1.back( ) = m1head;
          p1 = MergeKmerPaths( m1 );    }
     if (have_back_gap2)
     {    back_gap2 = m2.back( ).LastSegment( );
          KmerPath m2head;
          for ( int u = 0; u < m2.back( ).NSegments( ) - 1; u++ )
               m2head.AddSegmentNoConcatenate( m2.back( ).Segment(u) );
          m2.back( ) = m2head;
          p2 = MergeKmerPaths( m2 );    }

     // If p1 and/or p2 have adjacent and mergeable KmerPathIntervals, they will not
     // be handled correctly by ImproperMerge.  Therefore we contract mergeable
     // KmerPathIntervals in p1 and p2, and translate locations to there, reversing
     // the process later.

     KmerPath p1c, p2c;
     p1c.Reserve( p1.NSegments( ) );
     p2c.Reserve( p2.NSegments( ) );
     int w1c_index = w1.GetIndex( ), w2c_index = w2.GetIndex( );
     int w1c_loc = w1.GetLoc( ), w2c_loc = w2.GetLoc( );
     for ( int u = 0; u < p1.NSegments( ); u++ )
     {    int segs_before = p1c.NSegments( );
          int kmers_before = 0;
          if ( segs_before > 0 )
               kmers_before = p1c.Segment( segs_before - 1 ).Length( );
          p1c.AddSegment( p1.Segment(u) );
          int segs_after = p1c.NSegments( );
          if ( segs_before == segs_after )
          {    if ( u <= w1.GetIndex( ) ) --w1c_index;
               if ( u == w1.GetIndex( ) ) w1c_loc += kmers_before;    }    }
     for ( int u = 0; u < p2.NSegments( ); u++ )
     {    int segs_before = p2c.NSegments( );
          int kmers_before = 0;
          if ( segs_before > 0 )
               kmers_before = p2c.Segment( segs_before - 1 ).Length( );
          p2c.AddSegment( p2.Segment(u) );
          int segs_after = p2c.NSegments( );
          if ( segs_before == segs_after )
          {    if ( u <= w2.GetIndex( ) ) --w2c_index;
               if ( u == w2.GetIndex( ) ) w2c_loc += kmers_before;    }    }
     KmerPathLoc w1c( p1c, w1c_index, w1c_loc );
     KmerPathLoc w2c( p2c, w2c_index, w2c_loc );
     ForceAssertEq( w1c.GetKmer( ), w2c.GetKmer( ) );

     // Try to merge.
     vec<ImproperMerger> merges;
     {
       vec<ImproperMerger> merges0;
       ImproperMergePaths( p1c, p2c, w1c.GetIndex( ), w2c.GetIndex( ),
			   merges0, 1, ngv );
       for ( int u = 0; u < merges0.isize( ); u++ )
	 if ( ngv == NULL || ngv->MakeValid( merges0[u].merged ) )
	   merges.push_back( merges0[u] );

       if ( merges.empty( ) ) return 0;
     }
     M = merges[0].merged;

     // Translate positions on p1c and p2c back to p1 and p2, by counting kmers.

     KmerPathLoc left1c = merges[0].left_end1, right1c = merges[0].right_end1;
     KmerPathLoc left2c = merges[0].left_end2, right2c = merges[0].right_end2;
     int left1_count = 0, left2_count = 0, right1_count = 0, right2_count = 0;
     for ( int u = 0; u < left1c.GetIndex( ); u++ )
          left1_count += MinLen( p1c.Segment(u) );
     left1_count += left1c.GetLoc( );
     for ( int u = 0; u < left2c.GetIndex( ); u++ )
          left2_count += MinLen( p2c.Segment(u) );
     left2_count += left2c.GetLoc( );
     for ( int u = 0; u < right1c.GetIndex( ); u++ )
          right1_count += MinLen( p1c.Segment(u) );
     right1_count += right1c.GetLoc( );
     for ( int u = 0; u < right2c.GetIndex( ); u++ )
          right2_count += MinLen( p2c.Segment(u) );
     right2_count += right2c.GetLoc( );
     int left1_index = 0, left2_index = 0, right1_index = 0, right2_index = 0;
     int left1_loc = 0, left2_loc = 0, right1_loc = 0, right2_loc = 0;
     for ( int u = 0; u < p1.NSegments( ); u++ )
     {    if ( left1_count == MinLen( p1.Segment(u) ) )
          {    left1_loc = 0;
               left1_index = u+1;
               break;     }
          if ( left1_count < MinLen( p1.Segment(u) ) )
          {    left1_loc = left1_count;
               left1_index = u;
               break;    }
          left1_count -= MinLen( p1.Segment(u) );    }
     for ( int u = 0; u < p2.NSegments( ); u++ )
     {    if ( left2_count == MinLen( p2.Segment(u) ) )
          {    left2_loc = 0;
               left2_index = u+1;
               break;     }
          if ( left2_count < MinLen( p2.Segment(u) ) )
          {    left2_loc = left2_count;
               left2_index = u;
               break;    }
          left2_count -= MinLen( p2.Segment(u) );    }
     for ( int u = 0; u < p1.NSegments( ); u++ )
     {    if ( right1_count == MinLen( p1.Segment(u) ) )
          {    right1_loc = 0;
               right1_index = u+1;
               break;     }
          if ( right1_count < MinLen( p1.Segment(u) ) )
          {    right1_loc = right1_count;
               right1_index = u;
               break;    }
          right1_count -= MinLen( p1.Segment(u) );    }
     for ( int u = 0; u < p2.NSegments( ); u++ )
     {    if ( right2_count == MinLen( p2.Segment(u) ) )
          {    right2_loc = 0;
               right2_index = u+1;
               break;    }
          if ( right2_count < MinLen( p2.Segment(u) ) )
          {    right2_loc = right2_count;
               right2_index = u;
               break;    }
          right2_count -= MinLen( p2.Segment(u) );    }

     // Compute positions.

     KmerPathLoc left1( p1, left1_index, left1_loc );
     KmerPathLoc right1( p1, right1_index, right1_loc );
     KmerPathLoc left2( p2, left2_index, left2_loc );
     KmerPathLoc right2( p2, right2_index, right2_loc );
     ends.resize(2);
     ends[0].first = MultiKmerPathLoc( m1, left1 );
     ends[0].second = MultiKmerPathLoc( m1, right1 );
     ends[1].first = MultiKmerPathLoc( m2, left2 );
     ends[1].second = MultiKmerPathLoc( m2, right2 );
     ForceAssertEq( ends[0].first.Kmer(m1), ends[1].first.Kmer(m2) );
     ForceAssertEq( ends[0].second.Kmer(m1), ends[1].second.Kmer(m2) );
     i.resize( p.Path1( ).NVertices( ) ); j.resize( p.Path2( ).NVertices( ) );

     // If the calculate_ij flag was not set, skip the following chunk of code,
     // which is time-consuming and has no effect except on the output variables
     // i and j.
     if ( calculate_ij ) {
     KmerPath mc;
     for ( int x = 0; x < 2; x++ )
     {    vec<KmerPathSplitPoint>& ij = ( x == 0 ? i : j );
          for ( int u = 0; u < p.Path(x).NVertices( ); u++ )
          {    MultiKmerPathLoc uloc( u, 0, 0 );
               if ( uloc >= ends[x].first && uloc <= ends[x].second )
               {    KmerPathLoc l, onm;
                    if ( x == 0 ) onm = uloc.PosOnMerger(m1, p1);
                    else onm = uloc.PosOnMerger(m2, p2);
                    int index = onm.GetIndex( ), loc = onm.GetLoc( );
                    const KmerPath& p12 = ( x == 0 ? p1 : p2 );
		    mc.Clear( );
		    mc.Reserve( p12.NSegments( ) );
                    for ( int v = 0; v < p12.NSegments( ); v++ )
                    {    int segs_before = mc.NSegments( );
                         int kmers_before = 0;
                         if ( segs_before > 0 )
                              kmers_before = mc.Segment( segs_before - 1 ).Length( );
                         mc.AddSegment( p12.Segment(v) );
                         int segs_after = mc.NSegments( );
                         if ( segs_before == segs_after )
                         {    if ( v <= onm.GetIndex( ) ) --index;
                              if ( v == onm.GetIndex( ) )
                                   loc += kmers_before;    }    }
                    const KmerPath& p12c = ( x == 0 ? p1c : p2c );
                    KmerPathLoc onmc( p12c, index, loc );
                    if ( p12c.Segment(index).isGap( ) )
                    {    ForceAssertGe( index, 1 );
                         ForceAssertEq( loc, 0 );
                         ForceAssert( !p12c.Segment(index-1).isGap( ) );
                         KmerPathLoc onmcleft( p12c, index-1,
                              p12c.Segment(index-1).Length( ) - 1 );
                         KmerPathLoc lleft;
                         if ( x == 0 ) lleft = merges[0].PushLoc1(onmcleft);
                         else lleft = merges[0].PushLoc2(onmcleft);
                         ForceAssert( InRange(lleft) ); // XXX
                         Advance(lleft);
                         l = lleft;    }
                    else
                    {    if ( x == 0 ) l = merges[0].PushLoc1(onmc);
                         else l = merges[0].PushLoc2(onmc);
		    }
                    ij[u] = KmerPathSplitPoint( l.GetIndex( ), l.GetLoc( ) );    }
               const vecKmerPath& m12 = ( x == 0 ? m1 : m2 );
               if ( ends[x].second.AtRightEndOfEdge(m12)
                    && u == ends[x].second.Edge( ) + 1 )
               {    ij[u] = KmerPathSplitPoint(
                         M.NSegments( ), 0 );    }    }    }
     }

     // Restore front and back gaps if we had to remove them.  We try to extend the
     // alignment over the gaps, but we could try harder.

     MultiKmerPathLoc &e0a = ends[0].first, &e1a = ends[1].first;
     MultiKmerPathLoc &e0b = ends[0].second, &e1b = ends[1].second;
     if ( have_front_gap1 || have_front_gap2 )
     {    Bool equal_front
               = have_front_gap1 && have_front_gap2 && front_gap1 == front_gap2;
          if (equal_front)
          {
               // Check to see if the alignment goes right up to the gap on both
               // paths.

               if ( e0a.LeftEnd( ) && e1a.LeftEnd( ) )
               {
                    for ( int u = 0; u < 2; u++ )
                    {    if ( ends[u].second.Edge( ) == 0 )
                              ends[u].second.AddToSegmentOnEdge(1);    }
                    KmerPath Mnew;
                    Mnew.AddSegment(front_gap1);
                    for ( int u = 0; u < M.NSegments( ); u++ )
                         Mnew.AddSegmentNoConcatenate( M.Segment(u) );
                    M = Mnew;
                    for ( int u = 0; u < i.isize( ); u++ )
                    {    if ( i[u].IntervalsToLeft( ) > 0
                              || i[u].KmersInPartialToLeft( ) > 0 )
                         {    i[u].SetIntervalsToLeft(
                                   i[u].IntervalsToLeft( ) + 1 );    }    }
                    for ( int u = 0; u < j.isize( ); u++ )
                    {    if ( j[u].IntervalsToLeft( ) > 0
                              || j[u].KmersInPartialToLeft( ) > 0 )
                         {    j[u].SetIntervalsToLeft(
                                   j[u].IntervalsToLeft( ) + 1 );    }    }    }
               else
               {    if ( e0a.Edge( ) == 0 ) e0a.AddToSegmentOnEdge(1);
                    if ( e1a.Edge( ) == 0 ) e1a.AddToSegmentOnEdge(1);
                    if ( e0b.Edge( ) == 0 ) e0b.AddToSegmentOnEdge(1);
                    if ( e1b.Edge( ) == 0 ) e1b.AddToSegmentOnEdge(1);    }    }
          if ( !equal_front && have_front_gap1 )
          {    if ( e0a.Edge( ) == 0 ) e0a.AddToSegmentOnEdge(1);
               if ( e0b.Edge( ) == 0 ) e0b.AddToSegmentOnEdge(1);    }
          if ( !equal_front && have_front_gap2 )
          {    if ( e1a.Edge( ) == 0 ) e1a.AddToSegmentOnEdge(1);
               if ( e1b.Edge( ) == 0 ) e1b.AddToSegmentOnEdge(1);    }    }
     if ( have_back_gap1 || have_back_gap2 )
     {    Bool equal_back
               = have_back_gap1 && have_back_gap2 && back_gap1 == back_gap2;
          if (equal_back)
          {
               // Check to see if the alignment goes right up to the gap on both
               // paths.

               if ( static_cast<size_t>(e0b.Edge( )) == m1.size( ) - 1
                    && e0b.SegmentOnEdge( ) == m1.back( ).NSegments( ) - 1
                    && e0b.KmerOnSegment( )
                         == m1.back( ).LastSegment( ).Length( ) - 1
                    && static_cast<size_t>(e1b.Edge( )) == m2.size( ) - 1
                    && e1b.SegmentOnEdge( ) == m2.back( ).NSegments( ) - 1
                    && e1b.KmerOnSegment( )
                         == m2.back( ).LastSegment( ).Length( ) - 1 )
               {    M.AddSegment(back_gap1);
                    int ne1 = p.Path1( ).NEdges( ) - 1;
                    int ne2 = p.Path2( ).NEdges( ) - 1;
                    e0b.SetEdge(ne1);
                    e0b.SetSegmentOnEdge(
                         p.Path1( ).EdgeObject(ne1).NSegments( ) - 1 );
                    e0b.SetKmerOnSegment( p.Path1( ).EdgeObject(ne1).LastSegment( )
                         .Length( ) - 1 );
                    e1b.SetEdge(ne2);
                    e1b.SetSegmentOnEdge(
                         p.Path2( ).EdgeObject(ne2).NSegments( ) - 1 );
                    e1b.SetKmerOnSegment( p.Path2( ).EdgeObject(ne2).LastSegment( )
                         .Length( ) - 1 );
                    i.back( ) = j.back( ) = KmerPathSplitPoint( M.NSegments( ), 0 );
	       }    }    }
     return merges.isize( );
}







// InternalMergeImpl helper functions.

KmerPath Between( const KmerPath& p, const KmerPathSplitPoint& left,
     const KmerPathSplitPoint& right )
{    KmerPath answer;
     ForceAssert( left < right );                               
     ForceAssertLe( right.IntervalsToLeft( ), p.NSegments( ) );
     if ( left.IntervalsToLeft( ) == right.IntervalsToLeft( ) )
     {    const KmerPathInterval& I = p.Segment( left.IntervalsToLeft( ) );
          answer.AddSegment( I.Start( ) + left.KmersInPartialToLeft( ),
               I.Start( ) + right.KmersInPartialToLeft( ) - 1 );    }
     else
     {    const KmerPathInterval& Ileft = p.Segment( left.IntervalsToLeft( ) );
          if ( left.KmersInPartialToLeft( ) == 0 ) answer.AddSegment(Ileft);
          else
          {    answer.AddSegment( Ileft.Start( ) + left.KmersInPartialToLeft( ), 
                    Ileft.Stop( ) );    }
          for ( int u = left.IntervalsToLeft( ) + 1;  
               u < right.IntervalsToLeft( ); u++ )
          {    answer.AddSegment( p.Segment(u) );    }
          if ( right.KmersInPartialToLeft( ) > 0 )
          {    const KmerPathInterval& Iright 
                    = p.Segment( right.IntervalsToLeft( ) );
               answer.AddSegment( Iright.Start( ), Iright.Start( ) 
                    + right.KmersInPartialToLeft( ) - 1 );    }    }
     return answer;    }

void PrintLengths( const KmerPath& p )
{    for ( int i = 0; i < p.NSegments( ); i++ )
     {    if ( i > 0 ) cout << " ";
          cout << p.Segment(i).Length( );    }
     cout << "\n";    }

void PrintLengths( const EmbeddedSubPath<KmerPath> & p )
{    for ( int i = 0; i < p.NEdges( ); i++ )
     {    if ( i > 0 ) cout << " / ";
          for ( int j = 0; j < p.EdgeObject(i).NSegments( ); j++ )
          {    if ( j > 0 ) cout << " ";
               cout << p.EdgeObject(i).Segment(j).Length( );    }    }
     cout << "\n";    }

// SplitEdgeInPointedSubPath: given an EmbeddedSubPath u = (a,e) of a HyperKmerPath 
// 
//    e0       em-1
// a0 --> ...  --> am,
//
// and a kmer at position x on the EmbeddedSubPath (comprising a "PointedSubPath"), 
// an edge ej (either e0 or em-1), and a splitting of that edge in two (given by two 
// KmerPaths s1, s2 whose concatenation is the edge object ej), do the following:
// - edit the HyperKmerPath by splitting the edge;
// - update the data (u,x) accordingly.
// We assume that if the edge split point is to the left of x, then j = 0, and if
// the edge split point is to the right of x, then j = m-1 (and split_to_left_of_x
// must be specified, so giving j is redundant).  Moreover, in the first case, we
// edit u so that it starts at the edge split point, and in the second case, we
// edit u so that it stops at the edge split point.

void SplitEdgeInPointedSubPath( HyperKmerPath& h, EmbeddedSubPath<KmerPath>& u,
     MultiKmerPathLoc& x, int j, const KmerPath& s1, const KmerPath& s2,
     Bool split_to_left_of_x )
{    int segments_in_split_edge = u.EdgeObject(j).NSegments( );
     int next_to_last_vertex = u.Vertex( u.NVertices( ) - 2 );
     int new_vertex = h.N( );
     int new_edge1 = h.EdgeObjectCount( ), new_edge2 = h.EdgeObjectCount( ) + 1;
     h.SplitEdge( u.Vertex(j), u.EdgeObjectFromIndex(j), s1, s2 );
     int new_edge2_from_new_vertex 
          = h.EdgeObjectIndexToFromIndex( new_vertex, new_edge2 );
     int nsegs1 = s1.NSegments( ), nsegs2 = s2.NSegments( );
     Bool split_at_segment_boundary = ( nsegs1 + nsegs2 == segments_in_split_edge );
     if (split_to_left_of_x)
     {    ForceAssertEq( j, 0 );
          u.SetVertex( 0, new_vertex );
          u.SetEdge( 0, new_edge2_from_new_vertex );
          if ( x.Edge( ) == 0 )
          {    if ( !split_at_segment_boundary )
               {    x.SetSegmentOnEdge( x.SegmentOnEdge( ) - ( nsegs1 - 1 ) );
                    if ( x.SegmentOnEdge( ) == 0 )
                    {    x.SetKmerOnSegment( 
                              x.KmerOnSegment( ) 
                                   - s1.Segment( nsegs1 - 1 ).Length( ) );    }    }
               else x.SetSegmentOnEdge( x.SegmentOnEdge( ) - nsegs1 );    }    }
     else
     {    ForceAssertEq( j, u.NEdges( ) - 1 );
          int new_edge1_from_next_to_last_vertex
               = h.EdgeObjectIndexToFromIndex( next_to_last_vertex, new_edge1 );
          u.SetVertex( u.NVertices( ) - 1, new_vertex );
          u.SetEdge( u.NEdges( ) - 1, new_edge1_from_next_to_last_vertex );    }    }

// We define a class which can be stuffed into a priority queue and used to track
// extensions as we iteratively extend them through the graph.

class AlignExtend {
     
     public:

     AlignExtend( ) { }
     AlignExtend( const PointedSubPathPair& p, 
          const vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> >& ends,
          int align_length )
          : p_(p), ends_(ends), align_length_(align_length), terminal_(False) { }

     const PointedSubPathPair& PathPair( ) const { return p_; }

     const EmbeddedSubPath<KmerPath>& Path( int i ) const { return p_.Path(i); }
     const EmbeddedSubPath<KmerPath>& Path1( ) const { return p_.Path1( ); }
     const EmbeddedSubPath<KmerPath>& Path2( ) const { return p_.Path2( ); }

     const MultiKmerPathLoc& Loc1( ) const { return p_.Loc1( ); }
     const MultiKmerPathLoc& Loc2( ) const { return p_.Loc2( ); }

     const vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> >& Ends( ) const { return ends_; }
     const pair<MultiKmerPathLoc,MultiKmerPathLoc>& Ends(int i) const { return ends_[i]; }

     int AlignLength( ) const { return align_length_; }

     Bool Terminal( ) const { return terminal_; }

     void SetTerminal( const Bool term = True ) { terminal_ = term; }

     // Priority function for extensions: 
     // * alignments that can't be extended are lower priority that ones that can be.
     // * alignments that are longer are lower priority than ones that are shorter.
     // * alignments are then ordered by the vertices in their paths, so that alignments 
     //   of equal termination state, equal length, and equal start and vertices on both 
     //   components are considered equivalent.
     // This results in a breadth-first search where equivalent length paths that
     // terminate at the same set of vertices collapse into a single result (regardless
     // of the content of the path).  This compromise makes searching through diploid
     // components much faster.
     friend bool operator<( const AlignExtend& ae1, const AlignExtend& ae2 )
     {
       // Return false if ae1 is lower priority than ae2.
       if ( ae1.Terminal( ) && !ae2.Terminal( ) ) return false;
       if ( !ae1.Terminal( ) && ae2.Terminal( ) ) return true;
       if ( ae1.AlignLength( ) > ae2.AlignLength( ) ) return false;
       if ( ae1.AlignLength( ) < ae2.AlignLength( ) ) return true;
       for ( int p = 0; p < 2; ++p ) {
         if ( ae1.Path(p).FirstVertex() < ae2.Path(p).FirstVertex() ) return false;
         if ( ae1.Path(p).FirstVertex() > ae2.Path(p).FirstVertex() ) return true;
         if ( ae1.Path(p).LastVertex() < ae2.Path(p).LastVertex() ) return false;
         if ( ae1.Path(p).LastVertex() > ae2.Path(p).LastVertex() ) return true;
       }
       return false;
     }

     private:

     PointedSubPathPair p_;
     vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> > ends_;
     int align_length_;
     Bool terminal_;

};



// InternalMergeLeftExtend and InternalMergeRightExtend are helpers for 
// InternalMerge.  Given a PointedSubPathPair p and ends of an alignment of it, it 
// tries to enlarge the alignment by extending to the left [resp. right].  

void InternalMergeLeftExtend( const HyperKmerPath& h, set<AlignExtend>& Q,
			      const int max_Q_size,
     const NegativeGapValidator* ngv )
{    KmerPath M;
     vec<KmerPathSplitPoint> il, jl;
     vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> > endsplus;
     while( !Q.empty( ) && (int) Q.size( ) < max_Q_size ) {
       AlignExtend E = *( Q.begin( ) );
       if ( E.Terminal( ) ) break;
          Q.erase( Q.begin( ) );

          // If alignment does not go to the left end on either sequence, there is
          // nothing to do.

          if ( !E.Ends(0).first.LeftEnd( ) && !E.Ends(1).first.LeftEnd( ) )
          {    E.SetTerminal( );
               Q.insert(E);
               continue;    }
     
          // Handle case where alignment goes to the left end on one sequence
          // but not on the other.
          
          Bool found_extension = False;
          Bool not_to_end = 
               !E.Ends(0).first.LeftEnd( ) || !E.Ends(1).first.LeftEnd( );
          for ( int f = 0; f < 2; f++ )
          {    int g = 1 - f;
               if ( E.Ends(f).first.LeftEnd( ) && !E.Ends(g).first.LeftEnd( ) )
               {    int w = E.Path(f).Vertex(0);
                    for ( int i = 0; i < h.To(w).isize( ); i++ )
                    {    int v = h.To(w)[i];
                         const int edge = h.EdgeObjectIndexByIndexTo( w, i );
                         int ei = h.InputToOutputFrom( w, i );

			 // If the new edge already exists in either sequence,
			 // ignore this extension, rather than creating a loop or
			 // a crossing of paths.
			 EmbeddedSubPath<KmerPath> path1 = E.Path1( ), path2 = E.Path2( );
			 if ( path1.Contains( edge ) || path2.Contains( edge ) ) continue;

                         // Modify path by prepending e to the one sequence.
                         MultiKmerPathLoc loc1 = E.Loc1( ), loc2 = E.Loc2( );
                         if ( f == 0 )
                         {    path1.Prepend( v, ei );
                              loc1.SetEdge( loc1.Edge( ) + 1 );    }
                         else
                         {    path2.Prepend( v, ei );
                              loc2.SetEdge( loc2.Edge( ) + 1 );    }
			 
                         PointedSubPathPair pplus( path1, path2, loc1, loc2 );

                         // Align.

                         int naligns = AlignSubpaths( pplus, endsplus, M, il, jl, ngv );
                         if ( naligns == 0 ) continue;

                         // If the extension did not use the added edge, forget it.

                         if ( endsplus[f].first.Edge( ) != 0 ) continue;

                         // Save extension.

                         found_extension = True;
			 Q.insert( AlignExtend( pplus, endsplus, M.KmerCount( ) ) );
		    }
	       }
	  }
          if (not_to_end)
          {    if ( !found_extension )
               {    E.SetTerminal( );
                    Q.insert(E);    }
               continue;    }
     
          // Handle case where alignment goes to the left end on both sequences.

          int w1 = E.Path1( ).Vertex(0), w2 = E.Path2( ).Vertex(0);
	  int w1_size = h.To(w1).isize( );
          for ( int i1 = 0; i1 < w1_size; i1++ )
          {    int v1 = h.To(w1)[i1];
               int ei1 = h.InputToOutputFrom( w1, i1 );
	       int w2_size = h.To(w2).isize( );
               for ( int i2 = 0; i2 < w2_size; i2++ )
               {    int v2 = h.To(w2)[i2];
		    int edge1 = h.EdgeObjectIndexByIndexTo( w1, i1 );
		    int edge2 = h.EdgeObjectIndexByIndexTo( w2, i2 );
                    if ( edge1 == edge2 ) continue;

		    // If the new edges already exist in either sequence,
		    // ignore this extension, rather than creating a loop or
		    // a crossing of paths.
                    EmbeddedSubPath<KmerPath> path1 = E.Path1( ), path2 = E.Path2( );
		    if ( path1.Contains( edge1 ) || path2.Contains( edge1 ) ||
			 path1.Contains( edge2 ) || path2.Contains( edge2 ) ) continue;

                    // Modify p by prepending ei1 and ei2 to the sequences.
                    int ei2 = h.InputToOutputFrom( w2, i2 );
                    path1.Prepend( v1, ei1 ), path2.Prepend( v2, ei2 );
		    
                    MultiKmerPathLoc loc1 = E.Loc1( ), loc2 = E.Loc2( );
                    loc1.SetEdge( loc1.Edge( ) + 1 );
                    loc2.SetEdge( loc2.Edge( ) + 1 );
                    PointedSubPathPair pplus( path1, path2, loc1, loc2 );

                    // Align.

                    int naligns = AlignSubpaths( pplus, endsplus, M, il, jl, ngv );
                    if ( naligns == 0 ) continue;

                    // If the extension did not use the added edges, forget it.

                    if ( endsplus[0].first.Edge( ) != 0 ) continue;

                    // Save extension.

                    found_extension = True;
		    Q.insert( AlignExtend( pplus, endsplus, M.KmerCount( ) ) );
	       }
	  }
          if ( !found_extension )
          {    E.SetTerminal( );
               Q.insert(E);    }    }    }

void InternalMergeRightExtend( const HyperKmerPath& h, set<AlignExtend>& Q,
			       const int max_Q_size,
     const NegativeGapValidator* ngv )
{    KmerPath M;
     vec<KmerPathSplitPoint> il, jl;
     vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> > endsplus;
     while( !Q.empty( ) && (int) Q.size( ) < max_Q_size ) {
       AlignExtend E = *( Q.begin( ) );
       if ( E.Terminal( ) ) break;
          Q.erase( Q.begin( ) );

          // If alignment does not go to the right end on either sequence, there is
          // nothing to do.

          if ( !E.Ends(0).second.RightEnd( E.Path1( ) ) 
               && !E.Ends(1).second.RightEnd( E.Path2( ) ) )
          {    E.SetTerminal( );
               Q.insert(E);
               continue;    }

          // Handle case where alignment goes to the right end on one sequence
          // but not on the other.

          Bool found_extension = False;
          Bool not_to_end = !E.Ends(0).second.RightEnd( E.Path1( ) ) 
               || !E.Ends(1).second.RightEnd( E.Path2( ) );
          for ( int f = 0; f < 2; f++ )
          {    int g = 1 - f;
               if ( E.Ends(f).second.RightEnd( E.Path(f) ) 
                    && !E.Ends(g).second.RightEnd( E.Path(g) ) )
               {    int v = E.Path(f).Vertex( E.Path(f).NVertices( ) - 1 );
		    int v_size = h.From(v).isize( );
                    for ( int i = 0; i < v_size; i++ )
                    {    int w = h.From(v)[i];
                         const int edge = h.EdgeObjectIndexByIndexFrom( v, i );
                         int ei = i;

			 // If the new edge already exists in either sequence,
			 // ignore this extension, rather than creating a loop or
			 // a crossing of paths.
			 EmbeddedSubPath<KmerPath> path1 = E.Path1( ), path2 = E.Path2( );
			 if ( path1.Contains( edge ) || path2.Contains( edge ) ) continue;

                         // Modify p by appending e to the one sequence.
                         if ( f == 0 ) path1.Append( w, ei );
                         else path2.Append( w, ei );

                         MultiKmerPathLoc loc1 = E.Loc1( ), loc2 = E.Loc2( );
                         PointedSubPathPair pplus( path1, path2, loc1, loc2 );

                         // Align.

                         int naligns = AlignSubpaths( pplus, endsplus, M, il, jl, ngv );
                         if ( naligns == 0 ) continue;

                         // If the extension did not use the added edge, forget it.

                         int last_edge = E.Path(f).NEdges( );
                         if ( endsplus[f].second.Edge( ) != last_edge ) continue;

                         // Save extension.

                         found_extension = True;
			 Q.insert( AlignExtend( pplus, endsplus, M.KmerCount( ) ) );
		    }
	       }
	  }
          if (not_to_end)
          {    if ( !found_extension )
               {    E.SetTerminal( );
                    Q.insert(E);    }
               continue;    }
     
          // Handle case where alignment goes to the right end on both sequences.

          int v1 = E.Path1( ).Vertex( E.Path1( ).NVertices( ) - 1 ); 
          int v2 = E.Path2( ).Vertex( E.Path2( ).NVertices( ) - 1 );
	  int v1_size = h.From(v1).isize( );
          for ( int i1 = 0; i1 < v1_size; i1++ )
          {    int w1 = h.From(v1)[i1];
               int ei1 = i1;
	       int v2_size = h.From(v2).isize( );
               for ( int i2 = 0; i2 < v2_size; i2++ )
               {    int w2 = h.From(v2)[i2];
		    int edge1 = h.EdgeObjectIndexByIndexFrom( v1, i1 );
		    int edge2 = h.EdgeObjectIndexByIndexFrom( v2, i2 );
                    if ( edge1 == edge2 ) continue;

		    // If the new edges already exist in either sequence,
		    // ignore this extension, rather than creating a loop or
		    // a crossing of paths.
                    EmbeddedSubPath<KmerPath> path1 = E.Path1( ), path2 = E.Path2( );
		    if ( path1.Contains( edge1 ) || path2.Contains( edge1 ) ||
			 path1.Contains( edge2 ) || path2.Contains( edge2 ) ) continue;

                    // Modify p by appending ei1 and ei2 to the sequences.
                    int ei2 = i2;
                    path1.Append( w1, ei1 ), path2.Append( w2, ei2 );
		    
                    MultiKmerPathLoc loc1 = E.Loc1( ), loc2 = E.Loc2( );
                    PointedSubPathPair pplus( path1, path2, loc1, loc2 );

                    // Align.

                    int naligns = AlignSubpaths( pplus, endsplus, M, il, jl, ngv );
                    if ( naligns == 0 ) continue;

                    // If the extension did not use the added edges, forget it.

                    if ( endsplus[0].second.Edge( ) != E.Path1( ).NEdges( ) ) 
                         continue;

                    // Save extension.

                    found_extension = True;
		    Q.insert( AlignExtend( pplus, endsplus, M.KmerCount( ) ) );
	       }
	  }
          if ( !found_extension )
          {    E.SetTerminal( );
               Q.insert(E);    }    }    }



void InternalMergeImpl(HyperKmerPath& hkp, 
                       const NegativeGapValidator* ngv, 
                       int min_overlap, 
                       int min_proper_overlap, 
                       int max_Q_size,
                       Bool seed_unique, 
                       Bool seed_unique_weak, 
                       const vec<tagged_rpint>& uniqdb)
{    
  hkp.TestValid( );
  vec<KmerPathCorrespondence> shares;

  
  TaskTimer timer;
  timer.StartWithTimeOut(100 * 60);  // time in seconds
  

  while(1) {
    Bool found_merge = False;

    // - Find connected components of HyperKmerPath.

    equiv_rel e( hkp.N( ) );
    for ( int i = 0; i < hkp.N( ); i++ )
      for ( int j = 0; j < hkp.From(i).isize( ); j++ )
        e.Join( i, hkp.From(i)[j] );
    vec<int> reps;
    e.OrbitReps(reps);
    int nreps = reps.size( );
    vec< vec<int> > components0(nreps);
    for ( int i = 0; i < nreps; i++ )
      e.Orbit( reps[i], components0[i] );
    vec<int> to_rep( hkp.N( ) );
    for ( int i = 0; i < nreps; i++ ) {
      for ( int j = 0; j < components0[i].isize( ); j++ )
        to_rep[ components0[i][j] ] = i;    
    }

    // - If HyperKmerPath has more than one component, append reversed copy of 
    //   HyperKmerPath.

    int nverts = hkp.N( );
    int nedges = hkp.EdgeObjectCount( );
    if ( nreps > 1 )
    {    
      HyperKmerPath rev_copy(hkp);
      rev_copy.Reverse( );
      for ( int i = 0; i < nedges; i++ ) {
        KmerPath p = hkp.EdgeObject(i);
        p.Reverse( );
        rev_copy.SetEdgeObject( i, p );    
      }
      hkp.Append(rev_copy);
    }

    // - Set up structure to track connected components and correspondence
    //   with reverse complements.

    vec< vec<int> > components = components0;
    if ( nreps > 1 )
    {
      components.append(components0);
      for ( int i = nreps; i < 2 * nreps; i++ ) {
        for ( int j = 0; j < components[i].isize( ); j++ )
          components[i][j] += nverts;
      }
    }
    vec<int> touched( components.size( ), False );

    // - Find shared kmers, and remove those for which both copies lie in 
    //   the reversed copy of the original HyperKmerPath.  Also remove those
    //   which are between an edge and its reverse.  Also swap order in shares
    //   if needed so that reverse comes second.

    SharedKmers( hkp, shares, min_overlap );
    int count = 0;
    for ( int i = 0; i < shares.isize( ); i++ ) {
      if ( ( shares[i].From1( ) < nverts || shares[i].From2( ) < nverts )
           && shares[i].Id1( ) != shares[i].Id2( ) + nedges
           && shares[i].Id2( ) != shares[i].Id1( ) + nedges ) {
        if ( shares[i].From1( ) >= nverts ) shares[i].Swap( );
        if ( count != i ) shares[count] = shares[i];
        ++count;
      }
    }
    shares.resize(count);
    if (timer.TimedOut()) break;

    // - If seed_unique = True, remove shares whose perfect extension does not 
    //   contain unique sequence, except for cases where we can "zipper up".

    if ( seed_unique ) {
      vec<int> to_left, to_right;
      vec<int> edgelength;
      hkp.ToLeft(to_left), hkp.ToRight(to_right);
      edgelength.resize( hkp.EdgeObjectCount( ) );
      for ( int i = 0; i < edgelength.isize( ); i++ )
        edgelength[i] = hkp.EdgeLength(i);
      vec<Bool> to_remove( shares.size( ), False );
      for ( int i = 0; i < shares.isize( ); i++ ) {
        const KmerPathCorrespondence& s = shares[i];
        int e1 = s.Id1( ), e2 = s.Id2( );
        int v1 = to_left[e1], w1 = to_right[e1];
        int v2 = to_left[e2], w2 = to_right[e2];
        KmerPathLoc scanStart1( hkp.EdgeObject(e1), s.Pos1( ).Interval( ),
                                s.Pos1( ).PosOnInterval( ) );
        KmerPathLoc scanStart2( hkp.EdgeObject(e2), s.Pos2( ).Interval( ),
                                s.Pos2( ).PosOnInterval( ) );

        KmerPathLoc rightScan1( scanStart1 ), rightScan2( scanStart2 );
        ScanRightPerfectMatch( rightScan1, rightScan2 );
        int right_ext = KmersInInterval( scanStart1, rightScan1 ) - 1;

        KmerPathLoc leftScan1( scanStart1 ), leftScan2( scanStart2 );
        ScanLeftPerfectMatch( leftScan1, leftScan2 );
        int left_ext = KmersInInterval( leftScan1, scanStart1 ) - 1;

        int pos1 = s.Pos1( ).PosOnPath( ), pos2 = s.Pos2( ).PosOnPath( );
        int len1 = edgelength[e1], len2 = edgelength[e2];
        int left1 = pos1 - left_ext, right1 = pos1 + right_ext + 1;
        int left2 = pos2 - left_ext, right2 = pos2 + right_ext + 1;
        Bool zipper_up = False;
        if ( v1 == v2 && left1 == 0 && left2 == 0 ) zipper_up = True;
        if ( w1 == w2 && right1 == len1 && right2 == len2 ) zipper_up = True;
        if ( seed_unique && !zipper_up ) {
          if ( !seed_unique_weak || ( !hkp.Source(v1) && !hkp.Source(v2)
                                      && !hkp.Sink(w1) && !hkp.Sink(w2) ) ) {
            Bool uniq = False;
            KmerPath p;
            hkp.EdgeObject(e1).CopySubpath( leftScan1, rightScan1, p );
            vec<longlong> places;
            for ( int j = 0; j < p.NSegments( ); j++ ) {
              Contains( uniqdb, p.Segment(j), places );
              if ( places.nonempty( ) ) {
                uniq = True;
                break;
              }
            }
            if ( !uniq ) to_remove[i] = True;
          }
        }
      }
      EraseIf( shares, to_remove );
    }
    if (timer.TimedOut()) break;

    // - Note components corresponding to shares.

    vec<int> comp1( shares.size( ) ), comp2( shares.size( ) );
    for ( int i = 0; i < shares.isize( ); i++ ) {
      int v1 = shares[i].From1( ), v2 = shares[i].From2( );
      if ( v1 < nverts ) comp1[i] = to_rep[v1];
      else               comp1[i] = nreps + to_rep[ v1 - nverts ];
      if ( v2 < nverts ) comp2[i] = to_rep[v2];
      else               comp2[i] = nreps + to_rep[ v2 - nverts ];
    }
    //cout << "\t\tN shares: " << shares.size( ) << endl;
    if (timer.TimedOut()) break;
       
    // - Go through the shares.

    for ( int u = 0; u < shares.isize( ); u++ ) {
      const KmerPathCorrespondence& c = shares[u];
	    
      // -- Check share to see if it refers to an edge that is no longer
      //    active, as a result of edits to the graph.
	    
      if ( c.FromIndex1( ) >= hkp.From( c.From1( ) ).isize( )
           || hkp.EdgeObjectIndexByIndexFrom( c.From1( ), c.FromIndex1( ) ) 
           != c.Id1( ) )
        continue;
      if ( c.FromIndex2( ) >= hkp.From( c.From2( ) ).isize( )
           || hkp.EdgeObjectIndexByIndexFrom( c.From2( ), c.FromIndex2( ) ) 
           != c.Id2( ) )
        continue;
      if ( hkp.From( c.From1( ) )[ c.FromIndex1( ) ] != c.To1( ) ) continue;
      if ( hkp.From( c.From2( ) )[ c.FromIndex2( ) ] != c.To2( ) ) continue;
	    
      // -- Check share to see if it refers to the same edge (twice).
	    
      if ( c.Id1( ) == c.Id2( ) ) continue;
	    
      // -- Don't consider alignments of a component to its reverse.
	    
      if ( Abs( comp1[u] - comp2[u] ) == nreps ) continue;
	    
      // -- Once a component has been touched, we can't use its reverse.
	    
      if ( nreps > 1 ) {
        int r1, r2;
        if ( comp1[u] < reps.isize( ) ) r1 = comp1[u] + reps.isize( );
        else r1 = comp1[u] - reps.isize( );
        if ( comp2[u] < reps.isize( ) ) r2 = comp2[u] + reps.isize( );
        else r2 = comp2[u] - reps.isize( );
        if ( touched[r1] || touched[r2] ) continue;
      }
	    
      // -- Define PointedSubPairPath p.
	    
      vec<int> aa(2), bb(2);
      aa[0] = c.From1( ), aa[1] = c.To1( );
      bb[0] = c.From2( ), bb[1] = c.To2( );
      vec<int> e(1), f(1);
      e[0] = c.FromIndex1( );
      f[0] = c.FromIndex2( );
      PointedSubPathPair p( EmbeddedSubPath<KmerPath>( hkp, aa, e ),
                            EmbeddedSubPath<KmerPath>( hkp, bb, f ),
                            MultiKmerPathLoc( 0, c.Pos1( ) ), 
                            MultiKmerPathLoc( 0, c.Pos2( ) ) );
	    
      KmerPath M;
      vec< pair<MultiKmerPathLoc,MultiKmerPathLoc> > ends;
      vec<KmerPathSplitPoint> il, jl;
      int naligns = AlignSubpaths( p, ends, M, il, jl, ngv );
      if ( naligns == 0 ) continue;
      if (timer.TimedOut()) break;

      // -- Try to left extend and right extend.
            
      set<AlignExtend> Q;
      Q.insert( AlignExtend( p, ends, M.KmerCount( ) ) );
      InternalMergeLeftExtend( hkp, Q, max_Q_size, ngv ); 
       
      {
        // --- Clear all the terminal flags on the AlignExtends.
        set<AlignExtend> Qreset;
        set<AlignExtend>::iterator Q_iter;
        for ( Q_iter = Q.begin( ); Q_iter != Q.end( ); ++Q_iter ) {
          AlignExtend E = *Q_iter;
          E.SetTerminal( False );
          Qreset.insert( E );
        }
        Q = Qreset;
      }
      InternalMergeRightExtend( hkp, Q, max_Q_size, ngv );
	    
      // -- Go through each extension.
      while( !Q.empty( ) ) {
        //n_extensions++;
        //if ( n_extensions >= MAX_N_EXTENSIONS ) break;
	      
        AlignExtend E = *( Q.begin( ) );
        Q.erase( Q.begin( ) );
        p = E.PathPair( );
        ends = E.Ends( );
	      
        // --- Check to see if aligning region is too short.  There are 
        //     three criteria:
        //     1. Regardless of length, we accept any alignment if the 
        //        vertices agree on one end and the nonaligning parts on that 
        //        end have the same length (allowing for gap stretching).
        //     2. If the alignment is proper (i.e. goes from source to sink)
        //        and at least min_proper_overlap bases, we accept it.
        //     3. If the alignment extends for at least min_overlap bases, we 
        //        accept it.
	      
        Bool accept_short = False;
        if ( p.Path1( ).Vertex(0) == p.Path2( ).Vertex(0) ) {
          vec<int> before_min(2, 0), before_max(2, 0);
          for ( int i = 0; i < 2; i++ ) {
            const MultiKmerPathLoc& l = ends[i].first;
            const KmerPath& kp = p.Path(i).EdgeObject(0);
            int seg = l.SegmentOnEdge( ); 
            int kmer = l.KmerOnSegment( );
            for ( int j = 0; j < seg; j++ ) {
              before_min[i] += kp.MinLength(j);
              before_max[i] += kp.MaxLength(j);
            }
            before_min[i] += kmer;
            before_max[i] += kmer;
          }
          if ( IntervalOverlap( before_min[0], before_max[0] + 1,
                                before_min[1], before_max[1] + 1 ) > 0 )
            accept_short = True;
        }
        // --- No need to do the following chunk of code if we're already set the accept_short flag
        if ( !accept_short ) {
          if ( p.Path1( ).Vertex( p.Path1( ).NVertices( ) - 1 )
               == p.Path2( ).Vertex( p.Path2( ).NVertices( ) - 1 ) ) {
            vec<int> after_min(2, 0), after_max(2, 0);
            for ( int i = 0; i < 2; i++ ) {
              const MultiKmerPathLoc& r = ends[i].second;
              int last = p.Path(i).NEdges( ) - 1;
              const KmerPath& kp = p.Path(i).EdgeObject(last);
              int seg = r.SegmentOnEdge( ); 
              int kmer = r.KmerOnSegment( );
              after_min[i] += kp.Length(seg) - kmer - 1;
              after_max[i] += kp.Length(seg) - kmer - 1;
              for ( int r = seg + 1; r < kp.NSegments( ); r++ ) {
                after_min[i] += kp.MinLength(r);
                after_max[i] += kp.MaxLength(r);
              }
            }
            if ( IntervalOverlap( after_min[0], after_max[0] + 1,
                                  after_min[1], after_max[1] + 1 ) > 0 )
              accept_short = True;
          }
        }

        if ( !accept_short ) {
          if ( E.AlignLength( ) < min_proper_overlap ) continue;
		
          Bool left_proper =
            ( ends[0].first.LeftEnd( ) &&
              hkp.Source( p.Path(0).FirstVertex() ) ) ||
            ( ends[1].first.LeftEnd( ) &&
              hkp.Source( p.Path(1).FirstVertex() ) );
          Bool right_proper =
            ( ends[0].second.RightEnd( p.Path(0) ) &&
              hkp.Sink( p.Path(0).LastVertex( ) ) ) ||
            ( ends[1].second.RightEnd( p.Path(1) ) &&
              hkp.Sink( p.Path(1).LastVertex( ) ) );
		
          Bool proper = left_proper && right_proper;
		
          // ---- Check for pseudoproper: all the vertices and edges of a
          //      given component are covered.
		
          Bool pseudoproper = False;
          if ( !proper ) {
            for ( int j = 0; j < 2; j++ ) {
              const EmbeddedSubPath<KmerPath>& P = p.Path(j);
              if ( ends[j].first.LeftEnd( ) &&
                   ends[j].second.RightEnd(P) ) {
                vec<int> vertices;
                vec<int> edges;
                vertices.reserve( P.NVertices( ) );
                edges.reserve( P.NEdges( ) );
                for ( int u = 0; u < P.NVertices( ); u++ )
                  vertices.push_back( P.Vertex(u) );
                for ( int u = 0; u < P.NEdges( ); u++ ) {
                  edges.push_back( P.EdgeObjectIndex(u) );
                }
                UniqueSort(vertices), UniqueSort(edges);
                if ( hkp.IsComplete( vertices, edges ) ) {
                  pseudoproper = True;
                  break;
                }
              }
            }
          }
		
          if ( !proper && !pseudoproper ) {
            if ( E.AlignLength( ) < min_overlap ) 
              continue;    
          }
        }
        if (timer.TimedOut()) break;
        
        // --- Rebuild the alignment data once again;
        //     this time, for the first time, capture il and jl.
        naligns = AlignSubpaths( p, ends, M, il, jl, ngv, true );
        ForceAssertGt( naligns, 0 );
        if (timer.TimedOut()) break;

	      
        // --- Add vertices if ends are not at the ends of edges, and adjust
        //     p, il, jl accordingly.
	      
        for ( int i = 0; i < 2; i++ ) {
          const MultiKmerPathLoc& l = ends[i].first;
          if ( !l.LeftEnd( ) ) {
            const KmerPath& kp = p.Path(i).EdgeObject(0);
            int seg = l.SegmentOnEdge( ); 
            int kmer = l.KmerOnSegment( );
            KmerPath before, after;
            for ( int j = 0; j < seg; j++ )
              before.AddSegment( kp.Segment(j) );
            if ( kmer > 0 ) {
              before.AddSegment( kp.Segment(seg).Start( ),
                                 kp.Segment(seg).Start( ) + kmer - 1 );
            }
            if ( kmer == 0 ) after.AddSegment( kp.Segment(seg) );
            else {
              after.AddSegment( kp.Segment(seg).Start( ) + kmer,
                                kp.Segment(seg).Stop( ) );
            }
            for ( int r = seg + 1; r < kp.NSegments( ); r++ )
              after.AddSegment( kp.Segment(r) );
            SplitEdgeInPointedSubPath( hkp, p.PathMutable(i), 
                                       p.LocMutable(i), 0, before, after, True );    
            p.PathMutable(0).Repair( ), p.PathMutable(1).Repair( );
            if ( i == 0 ) il[0] = KmerPathSplitPoint( 0, 0 );
            if ( i == 1 ) jl[0] = KmerPathSplitPoint( 0, 0 );
            MultiKmerPathLoc& r = ends[i].second;
            if ( l.Edge( ) == r.Edge( ) ) {
              if ( l.SegmentOnEdge( ) == r.SegmentOnEdge( ) ) {
                r.SetKmerOnSegment( r.KmerOnSegment( ) 
                                    - l.KmerOnSegment( ) );
              }    
              r.SetSegmentOnEdge( r.SegmentOnEdge( ) 
                                  - l.SegmentOnEdge( ) );
            }
          }
          const MultiKmerPathLoc& r = ends[i].second;
          if ( !r.RightEnd( p.Path(i) ) ) {
            int last = p.Path(i).NEdges( ) - 1;
            const KmerPath& kp = p.Path(i).EdgeObject(last);
            int seg = r.SegmentOnEdge( ); 
            int kmer = r.KmerOnSegment( );
            KmerPathLoc right( kp, seg, kmer );
            KmerPath before, after;
            kp.CopyHead( right, before );
            if ( kmer < kp.Segment(seg).Length( ) - 1 )
              after.AddSegment( kp.Segment(seg).Start( ) + kmer + 1,
                                kp.Segment(seg).Stop( ) );
            for ( int r = seg + 1; r < kp.NSegments( ); r++ )
              after.AddSegment( kp.Segment(r) );
            SplitEdgeInPointedSubPath( hkp, p.PathMutable(i), p.LocMutable(i), 
                                       last, before, after, False );    
            p.PathMutable(0).Repair( ); 
            p.PathMutable(1).Repair( );    
            if ( i == 0 ) 
              il[ p.Path1( ).NVertices( ) - 1 ]
                = KmerPathSplitPoint( M.NSegments( ), 0 );
            if ( i == 1 ) 
              jl[ p.Path2( ).NVertices( ) - 1 ]
                = KmerPathSplitPoint( M.NSegments( ), 0 );
          }
        }
        if (timer.TimedOut()) break;
            
        // --- Glue.
        vec<KmerPathSplitPoint> splits = il;
        splits.append(jl);
        UniqueSort(splits);
        vec<int> EE( il.size( ) ), FF( jl.size( ) );
        for ( int v = 0; v < il.isize( ); v++ )
          EE[v] = BinPosition( splits, il[v] );
        for ( int v = 0; v < jl.isize( ); v++ )
          FF[v] = BinPosition( splits, jl[v] );
        vec<KmerPath> Msplit;
        for ( int v = 0; v < splits.isize( ) - 1; v++ )
          Msplit.push_back( Between( M, splits[v], splits[v+1] ) );
        digraphE<KmerPath> Mgraph( Msplit, digraphE<KmerPath>::EDGES_IN_LINE );
	      
        hkp.Glue( p.Path1( ), p.Path2( ), EE, FF, Mgraph );
        found_merge = True;
	      
        // --- Note which components have been touched.
	      
        touched[ comp1[u] ] = touched[ comp2[u] ] = True;
        break;
      }
      if (timer.TimedOut()) break;
    }
    if (timer.TimedOut()) break;
   
	  
    // - Remove superfluous components from original introduction of reverse
    //   complements.

    if ( nreps > 1 )
    {    
      for ( int i = 0; i < nreps; i++ ) {
        int del;
        ForceAssert( !touched[i] || !touched[i+nreps] );
        if ( !touched[i] && !touched[i+nreps] ) del = i + nreps;
        else if ( touched[i] ) del = i + nreps;
        else del = i;
        for ( int j = 0; j < components[del].isize( ); j++ ) {
          int v = components[del][j];
          hkp.DeleteEdgesAtVertex(v);
        }
      }
    }
    hkp.RemoveUnneededVertices( );    

    // - Remove identical paths between two given vertices.  Then remove dead
    //   edge objects.

    hkp.CompressEdgeObjects( );
    hkp.RemoveDuplicateEdges( );
    hkp.RemoveDeadEdgeObjects( );

    if ( !found_merge ) break;
  }


  hkp.TestValid( );
}

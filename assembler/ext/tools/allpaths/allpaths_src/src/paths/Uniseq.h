///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A 'uniseq' is a sequence of unibases, with defined overlap between them,
// usually K-1, but possibly smaller.
//
// A 'snark' is a kind of assembly in progress, that is a graph encoded in terms of 
// unibases, and which can contain gaps and ambiguities.
//
// More specifically, let a 'gapster' denote either a gap (represented by separation 
// and deviation), or a list of uniseq closures.  Then a snark is a directed graph 
// whose vertices are uniseq objects and whose edges are gapster objects.  A snark
// may have at most one edge between two vertices.  It is understood that within a 
// snark, the uniseqs in a closed gapster share one unibase with the connecting 
// uniseqs on either end, as in
//
//               (1,2,3) --{ (3,4,5,6), (3,4,5,4,5,6) }--> (6,7,8)
//               uniseq    ...... closed gapster ......    uniseq
//
// which could have come from an open gapster between two uniseqs:
//
//               (1,2,3) --------( 1000 +/- 100 )--------> (6,7,8).
//
// Note that you must call SetUnibases to set the static location of unibases, for
// the classes uniseq and snark.  Otherwise bad things will happen.

#ifndef UNISEQ_H
#define UNISEQ_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "graph/Digraph.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "paths/Sepdev.h"

class uniseq {

     public:

     uniseq( ) { }
     uniseq( const vec<int>& u, const vec<int>& overlap )
          : u_(u), overlap_(overlap)
     {    AssertEq( u.size( ), overlap.size( ) + 1 );    }

     Bool Dead( ) const { return u_.empty( ) && overlap_.empty( ); }
     Bool Alive( ) const { return !Dead( ); }
     void Kill( ) { u_.clear( ); overlap_.clear( ); }

     const vecbasevector& Unibases( ) const 
     {    Assert( unibases_ != 0 );
          return *unibases_;    }
     const basevector& Unibase( int k ) const 
     {    Assert( unibases_ != 0 );
          return Unibases( )[k];    }
     void SetUnibases( const vecbasevector& unibases ) { unibases_ = &unibases; }

     int N( ) const { return u_.size( ); }
     const vec<int>& U( ) const { return u_; }
     int U( int k ) const { return u_[k]; }
     const vec<int>& Over( ) const { return overlap_; }
     int Over( int k ) const { return overlap_[k]; }

     int Len( ) const; // return length in bases

     basevector Bases( ) const;

     // TrimLeft: remove n unibases from the left.
     // TrimRight: remove n unibases from the left.

     void TrimLeft( const int n );
     void TrimRight( const int n );
     
     void TrimEnds( const int n, const int m );

     void ReverseMe( const vec<int>& to_rc );
     uniseq Reverse( const vec<int>& to_rc ) const;

     friend ostream& operator<<( ostream& out, const uniseq& q );

     void Print( ostream& out, const int K ) const;

     friend void Print( ostream& out, const vec< vec<uniseq> >& ul, const int K );

     // Concatenate two uniseqs, assuming that they share a unibase in the middle
     // that should be deleted.  Ditto for three.

     friend uniseq Cat( const uniseq& s1, const uniseq& s2 );
     friend uniseq Cat( const uniseq& u1, const uniseq& u2, const uniseq& u3 );

     Bool Contains( const uniseq& x ) const;

     // Contains( x, p ).  You must have either p = 0 or p = -1.  If p = 0,
     // return True if *this contains x at its beginning.  If p = -1, return True
     // if *this contains x at its end.

     Bool Contains( const uniseq& x, int p ) const
     {    ForceAssert( p == 0 || p == -1 );
          if ( p == 0 )
               return U( ).Contains( x.U( ), 0 ) && Over( ).Contains( x.Over( ), 0 );
          else
          {    if ( !U( ).Contains( x.U( ), N( ) - x.N( ) ) ) return False;
               return Over( ).Contains( x.Over( ), N( ) - x.N( ) );    }    }

     friend vec<int> Common( const vec<uniseq>& v );

     friend Bool operator==( const uniseq& u1, const uniseq& u2 )
     {    return u1.U( ) == u2.U( ) && u1.Over( ) == u2.Over( );    }
     friend Bool operator!=( const uniseq& u1, const uniseq& u2 )
     {    return !( u1 == u2 );    }

     friend Bool operator<( const uniseq& u1, const uniseq& u2 )
     {    if ( u1.U( ) < u2.U( ) ) return True;
          if ( u1.U( ) > u2.U( ) ) return False;
          return u1.Over( ) < u2.Over( );    }

     private:

     vec<int> u_;       // unibases ids
     vec<int> overlap_; // overlap between successive unibases
     const static vecbasevector* unibases_;

};

class gapster {
  
     public:

     gapster( ) : open_(True), sep_(0), dev_(0) { }
     gapster( const int sep, const int dev ) : open_(True), sep_(sep), dev_(dev) { }

     gapster( const uniseq& u ) : open_(False), sep_(0), dev_(0) {
       closures_.push_back(u); 
     }

     gapster( vec<uniseq> c ) : open_(False), sep_(0), dev_(0) {
       closures_ = c; 
     } 

     Bool Open( ) const { return open_; }
     Bool Closed( ) const { return !open_; }

     void Close( const vec<uniseq>& closures )
     {    closures_ = closures;
          open_ = False;    }

     void RemoveSomeClosures( const vec<Bool>& to_delete )
     {    EraseIf( closures_, to_delete );    }

     int Sep( ) const
     {    Assert( Open( ) );
          return sep_;    }
     int Dev( ) const
     {    Assert( Open( ) );
          return dev_;    }

     int MinLen( const double dev_mult ) const;
     int MaxLen( const double dev_mult ) const;
     int MidLen( ) const 
     {    return int( round( double( MinLen(3.0) + MaxLen(3.0) ) / 2.0 ) );    }

     const vec<uniseq>& Closures( ) const { return closures_; }

     int ClosureCount( ) const // OK to call for open or closed
     {    return closures_.size( );    }

     const uniseq& Closure( int n ) const
     {    Assert( Closed( ) );
          return closures_[n];    }      

     friend Bool operator==( const gapster& g1, const gapster& g2 )
     {    if ( g1.open_ != g2.open_ ) return False;
          if ( g1.open_ ) return g1.sep_ == g2.sep_ && g1.dev_ == g2.dev_;
          else return g1.closures_ == g2.closures_;    }

     private:

     Bool open_;
     int sep_;
     int dev_;
     vec<uniseq> closures_;

};

class snark {

     public:

     snark( ) { }
     snark( const digraphE<gapster>& G, const vec<uniseq>& seq )
          : G_(G), seq_(seq) { }

     const digraphE<gapster>& G( ) const { return G_; }
     digraphE<gapster>& Gmutable( ) { return G_; }

     const vec<int>& From( int j ) const { return G_.From(j); }
     const vec<int>& To( int j ) const { return G_.To(j); }

     int VertN( ) const { return G( ).N( ); }
     int EdgeN( ) const { return G( ).EdgeObjectCount( ); }

     const uniseq& Vert( int k ) const { return seq_[k]; }
     uniseq& VertMutable( int k ) { return seq_[k]; }
     const gapster& Edge( int k ) const { return G( ).EdgeObject(k); }

     const vecbasevector& Unibases( ) const 
     {    Assert( unibases_ != 0 );
          return *unibases_;    }
     const basevector& Unibase( int u ) const 
     {    Assert( unibases_ != 0 );
          return Unibases( )[u];    }
     void SetUnibases( const vecbasevector& unibases ) { unibases_ = &unibases; }

     const vec<int>& ToRc( ) const
     {    Assert( to_rc_ != 0 );
          return *to_rc_;    }
     int ToRc( int u ) const
     {    Assert( to_rc_ != 0 );
          return ToRc( )[u];    }
     void SetToRc( const vec<int>& to_rc ) { to_rc_ = &to_rc; }

     void CloseGap( const int v1, const int v2, const vec<uniseq>& closures );

     // SwallowSimpleGaps finds and eliminates certain gaps having single closures.

     void SwallowSimpleGaps( );

     // Remove dead vertices:

     void BringOutTheDead( );

     // Simple inverted repeats.  Look for an assembly whose global structure
     // can be expressed as:
     // X ----> X'
     //   ---->
     //   <----
     //   <----.
     // In such cases, we can delete one of the arrows.  We delete one of the
     // *shorter* edges.  This might be generalized.

     void HandleSimpleInvertedRepeats( );

     void RemoveSubsumedStuff( );

     // EstimatedGenomeSize returns the sum of the lengths of the vertices and
     // the "mid" lengths of the edges.

     int64_t EstimatedGenomeSize( ) const;

     // Compute scaffolds

     void ComputeScaffolds( vec<superb>& superbs, vec<efasta>& econtigs, 
          digraphE<sepdev>& SG ) const;

     private:

     digraphE<gapster> G_;
     vec<uniseq> seq_;
     const static vecbasevector* unibases_;
     const static vec<int>* to_rc_;

};

// A placement_on describes the position of read on a snark.  When placed on a
// closed edge object, we do not track the closure that the read lies on, and thus
// a single placement_on may represent several possible positions.
//
// Note that 'id' is a vertex id if it's less than #vertices; otherwise
// id - #vertices is an edge id.

class placement_on {
     
     public:

     placement_on( ) { }
     placement_on( const int id, const int pos, const int apos, const Bool fw )
          : id(id), pos(pos), apos(apos), fw(fw) { }

     int id;    // object identifier (vertex if < #vertices, else edge, offset)
     int pos;   // start position on object, in bases
     int apos;  // object length - stop position on object
     Bool fw;   // placed forward?

     void Print( ostream& out, const snark& S, const vec<int>& to_left,
          const vec<int>& to_right );

};

// GetPairPlacements finds all placements of a read pair on a snark.  Note that it
// will not find placements where one or both ends dangle in open gaps or outside
// the snark entirely.  Note also that end placements in closed edges are tracked
// only by base position, without reference to the particular closure(s) on which
// the read lies.

void GetPairPlacements(

     // assembly:

     const snark& S, const vec<int>& to_right,
     const vec<int>& min_len, const vec<int>& max_len,

     // input pair:

     const int64_t pid,

     // reads and pairs:

     const vecbasevector& bases,
     const PairsManager& pairs,

     // read placements on unipaths -- read --> ( unipath, pos, fw? ):

     const vec< vec< triple<int,int,Bool> > >& placements_by_read,

     // edge_id --> ( u, pos, apos ):

     const vec< vec< triple<int,int,int> > >& u_pos_apos,

     // placements of the reads -- pos_on[0] provides placements of the first read
     // and pos_on[1] provided placements of the second read:

     vec< vec<placement_on> >& pos_on,

     // placements of the pairs -- the first and second entries are indices in
     // pos_on[0] and pos_on[1], respectively; the third entry is set if the second
     // read is fw rather than the first:

     vec< triple<int,int,Bool> >& placements );

void FindPartnersInGap(

     // vertices to look between

     const int x1, const int x2,

     // distance to scan on left and right

     const int flank,

     // assembly

     const snark& S,

     // reads and pairs

     const vecbasevector& bases, const PairsManager& pairs,

     // 1. read --> ( unipath, pos, fw? )
     // 2. unipath --> ( read_id, pos, fw? )

     const vec< vec< triple<int,int,Bool> > >& placements_by_read,
     const vec< vec< triple<int64_t,int,Bool> > >& placements_by_unipath,

     // ( unipath, start, stop )

     vec< triple<int,int,int> >& between );

#endif

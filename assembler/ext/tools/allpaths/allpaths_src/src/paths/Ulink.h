///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ULINK_H
#define ULINK_H

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"

class ulink_with_uids {

     public:

     ulink_with_uids( ) { }

     ulink_with_uids( const int u1, const int u2, const int sep, const int dev,
          const int start1, const int stop2, const int nlinks ) : u1(u1), u2(u2), 
          sep(sep), dev(dev), start1(start1), stop2(stop2), nlinks(nlinks) { }

     friend Bool operator<( const ulink_with_uids& x, const ulink_with_uids& y )
     {    if ( x.u1 < y.u1 ) return True;
          if ( x.u1 > y.u1 ) return False;
          if ( x.u2 < y.u2 ) return True;
          if ( x.u2 > y.u2 ) return False;
          if ( x.sep < y.sep ) return True;
          if ( x.sep > y.sep ) return False;
          if ( x.dev < y.dev ) return True;
          if ( x.dev > y.dev ) return False;
          if ( x.start1 < y.start1 ) return True;
          if ( x.start1 > y.start1 ) return False;
          if ( x.stop2 < y.stop2 ) return True;
          if ( x.stop2 > y.stop2 ) return False;
          if ( x.nlinks < y.nlinks ) return True;
          return False;    }

     int u1, u2;
     int sep, dev;
     int start1, stop2;
     int nlinks;

};

class ulink {

     public:

     ulink( ) { }
     ulink( const int sep, const int dev, const int start1, const int stop2,
          const int nlinks )
          : sep(sep), dev(dev), start1(start1), stop2(stop2), nlinks(nlinks) { }

     int sep, dev;
     int start1, stop2;
     int nlinks;

     int Sep( ) const { return sep; }
     int Dev( ) const { return dev; }
     int Var( ) const { return dev*dev; }

};

inline ulink CombineUlinks( const ulink& x, const ulink& y )
{    ulink z;
     vec<normal_distribution> N;
     N.push( x.sep, x.dev );
     N.push( y.sep, y.dev );
     normal_distribution n = CombineNormalDistributions(N);
     z.sep = int( round( n.mu_ ) );
     z.dev = int( round( n.sigma_ ) );
     z.nlinks = x.nlinks + y.nlinks;
     z.start1 = Min( x.start1, y.start1 );
     z.stop2 = Max( x.stop2, y.stop2 );
     return z;    }

// If we have edges a --> b and b --> a and the support for one direction
// is much stronger than the other, delete the other.

void DeleteBi( digraphE<ulink>& G, const int bi_mult, const Bool log );

// If we have edges a --> b and b --> a and c --> a and c --> b, and this
// evidence strongly supports one of a <--> b over the other, delete the other.

void DeleteBiAnother( digraphE<ulink>& G, const vecbasevector& unibases,
     const Bool log );

// Delete edges with nlinks < min_links EXCEPT in the situation where we have 
// e1: a --> b (nlinks >= min_links)
// e2: b --> c (nlinks < min_links)
// e3: a --> c (nlinks >= min_links).
// Some other exceptions should be added.

void DeleteNlinks( digraphE<ulink>& G, const int min_links, const Bool log );

// Delete transitive edges.  Transfer coverage from deleted edges.

void DeleteTrans( digraphE<ulink>& G, const vecbasevector& unibases,
     const int trans_depth, const int max_sep, const Bool log );

// Delete edges that are weakly supported.

void DeleteWeak( digraphE<ulink>& G, const vecbasevector& unibases,
     const int min_glue, const Bool log );

#endif

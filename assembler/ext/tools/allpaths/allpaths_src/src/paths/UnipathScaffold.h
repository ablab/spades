///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_SCAFFOLD_H
#define UNIPATH_SCAFFOLD_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Superb.h"
#include "graph/Digraph.h"

class linklet {

     public:

     linklet( ) { }
     linklet( const int sep, const int dev, const int nlinks,
          const int nlinks_indirect ) : sep(sep), dev(dev), nlinks(nlinks),
          nlinks_indirect(nlinks_indirect) {  }

     int sep;
     int dev;
     int nlinks;
     int nlinks_indirect;

     friend Bool operator<( const linklet& l1, const linklet& l2 )
     {    if ( l1.nlinks > l2.nlinks ) return True;
          if ( l1.nlinks < l2.nlinks ) return False;
          if ( l1.nlinks_indirect > l2.nlinks_indirect ) return True;
          if ( l1.nlinks_indirect < l2.nlinks_indirect ) return False;
          if ( l1.sep < l2.sep ) return True;
          if ( l1.sep > l2.sep ) return False;
          return l1.dev < l2.dev;    }

};

void UnipathScaffold(

     // input:

     const digraphE<linklet>& G,
     const vecbasevector& unibases,
     const int K,
     const vec<int>& to_rc,
     const vec<Bool>& exclude,
     const int min_kmers,

     // logging and debugging:

     const int FORCE_FIRST,
     const Bool SCAFFOLDING_VERBOSE,

     // output:

     vec<superb>& uscaffolds,
     vec<Bool>& circled );

void UnipathScaffoldAlt(

     // input:

     const digraphE<linklet>& G,
     const vecbasevector& unibases,
     const int K,
     const vec<int>& to_rc,
     const vec<Bool>& exclude,
     const int min_kmers,

     // logging and debugging:

     const int FORCE_FIRST,
     const Bool SCAFFOLDING_VERBOSE,

     // output:

     vec<superb>& uscaffolds,
     vec<Bool>& circled );

#endif

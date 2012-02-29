/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_WALK_H
#define SIMPLE_WALK_H

// SimpleWalkRight.  Try to close a given pair p using given reads, yielding
// closures.  The vector L contains the lengths, e.g. of p.Left(i).  If 
// max_opens > 0 and the number of open nodes in the graph comes to exceed
// max_opens, set fail = True and return no closures.

#include "CoreTools.h"
#include "paths/PairedPair.h"

void SimpleWalkRight( const pp_pair& p, const vec<pp_read>& reads, 
     const vec<int>& L, const double dmult, vec<pp_closure>& closures,
     const int max_opens, Bool create_closures_if_fail, Bool& fail, int verbosity,
     Bool depth_first = False, int max_nodes = 0 );

#endif

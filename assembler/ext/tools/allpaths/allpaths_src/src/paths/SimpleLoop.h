/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_LOOP_H
#define SIMPLE_LOOP_H

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"

class simple_loop {

     public:
     
     simple_loop( ) { }
     simple_loop( int u, int v, int w, int uv, int vv, int vw )
          : u(u), v(v), w(w), uv(uv), vv(vv), vw(vw) { }
          
     int u, v, w;     // vertices
     int uv, vv, vw;  // edges
     
};

void GetSimpleLoops( const HyperKmerPath& h, vec<simple_loop>& loops );

void DisambiguateSimpleLoops( HyperKmerPath& h, const vec<simple_loop>& loops,
     const vec< vec< pair<int,int> > >& sep_dev );

void DisambiguateSimpleLoops2( HyperKmerPath& h, HyperBasevector& hb,
     const vecbasevector& reads, const PairsManager& pairs,
     const vec<alignment_plus>& Aligns,
     const Bool verbose = False );


#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BreakBadScaffolds. 


#ifndef BREAK_BAD_SCAFFOLDS_H
#define BREAK_BAD_SCAFFOLDS_H

#include <omp.h>

#include "Fastavector.h"
#include "Charvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "math/NStatsTools.h"
#include "paths/Alignlet.h"


class blink_fw {

     public:

     int start1; // start on t1
     int stop1_low, stop1_high; // where read2 would land if it did
     int t2;
     int stop2;  // stop on t2

     friend Bool operator<( const blink_fw& b1, const blink_fw& b2 )
     {    if ( b1.t2 < b2.t2 ) return True;
          if ( b1.t2 > b2.t2 ) return False;
          if ( b1.start1 < b2.start1 ) return True;
          return False;    }

};

class blink_rc {

     public:

     int stop1; // stop on t1
     int start1_low, start1_high; // where read2 would land if it did
     int t2;
     int start2;  // start on t2

     friend Bool operator<( const blink_rc& b1, const blink_rc& b2 )
     {    if ( b1.t2 < b2.t2 ) return True;
          if ( b1.t2 > b2.t2 ) return False;
          if ( b1.stop1 < b2.stop1 ) return True;
          return False;    }

};

void ReportScaffoldsN50( const vec<superb> &supers, ostream &out );


void break_scaffolds( vec<superb> & scaffolds, vec<superb> & new_scaffolds,
		      PairsManager & pairs, vec<alignlet> & aligns0, vec<int> & aligns0_index, 
		      const int nreads, const int min_reach_away, vec<int> & trace_ids,
		      ofstream & log, vec<fastavector> & contigs, ofstream & sout, ofstream & cgout );

#endif

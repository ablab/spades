///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__PROCESS_GAP_H
#define PATHS__PROCESS_GAP_H

#include "Basevector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "paths/BigMapTools.h"
#include "paths/LongReadTools.h"
#include "paths/Ulink.h"
#include "paths/UnipathScaffold.h"
#include "paths/Uniseq.h"

void GetWalks( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc, 
     const vec< vec< pair<int,int> > >& nextsx, const vec<int>& use,
     vec< vec< pair<int,int> > >& walks1, int& bad );

void WalksToPatches( const vec< vec< pair<int,int> > >& walks1,
     const vecbasevector& unibases, vec<basevector>& patches );

void FilterWalksUsingJumps( const int u1, const int u2,
     const vec<basevector>& jbases_sorted, const vec<int64_t>& jbases_sorted_id,
     const PairsManager& jpairs, const vec< triple<int64_t,int,int> >& jaligns,
     const vecbasevector& unibases, vec<basevector>& patches,
     vec< vec< pair<int,int> > >& walks1 );

void FilterWalksUsingJumps2( const int edge_id, snark& S, const vec<int>& to_left,
     const vec<int>& to_right, const vecbasevector& jbases,    
     const PairsManager& jpairs,
     const vec< vec< triple<int,int,Bool> > >& placements_by_read,
     const vec< vec< triple<int64_t,int,Bool> > >& placements_by_unipath,
     const int verbosity );

void FilterWalksByIllegalCN1( const superb& s,
     vec<basevector>& patches, vec< vec< pair<int,int> > >& walks1 );

Bool Follow( const placementy& p1, const placementy& p2, const int K );

void Print( ostream& out, const placementy& p, const int len );

void ProcessGap(

     // definition of the gap

     const superb& s, const int u1, const int u2, const int sep, const int dev,

     // global structures

     const vecbasevector& unibases, const vec<int>& to_rc,
     const vec< vec< pair<int,int> > >& nextsx, 
     const vec< vec< pair<int,int> > >& nextsy, const int K, const vec<int>& use,

     // jump stuff to filter

     const vec<basevector>& jbases_sorted, const vec<int64_t>& jbases_sorted_id,
     const PairsManager& jpairs, const vec< triple<int64_t,int,int> >& jaligns,
     const vec< vec< triple<int,int,Bool> > >& placements_by_read,

     // for evaluation and logging

     ostream& out, const int verbosity, const Bool validate, const int LG, 
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs,
     Bool& have_lastp, placementy& lastp, int& gapcount, const digraphE<linklet>& G,

     // output

     Bool& found_patch, efasta& epatch, int& lastover,
     vec<int>& extras    // unibases1 used in every patch

);

#endif

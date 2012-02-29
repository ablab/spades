/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_paths_UnibaseUtils_h
#define __INCLUDE_paths_UnibaseUtils_h

#include "Basevector.h"
#include "CommonSemanticTypes.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/GetNexts.h"

/**
   Function: UnibaseInvolution

   For each <unibase>, identify its reverse complement.
 */
void UnibaseInvolution( const vecbasevector& unibases, vec< int >& toRc, nbases_t K );


/**
   Function: BuildUnibaseAdjacencyGraph

 */
void BuildUnibaseAdjacencyGraph( const vecbasevector& unibases, digraph& AG, nbases_t K );
#endif
// __INCLUDE_paths_UnibaseUtils_h

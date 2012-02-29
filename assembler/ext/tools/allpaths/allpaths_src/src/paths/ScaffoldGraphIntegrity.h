///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SCAFFOLD_GRAPH_INTEGRITY_H
#define SCAFFOLD_GRAPH_INTEGRITY_H

#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * ScaffoldGraphIntegrity
 *
 * Consistency tests for a scaffold graph.
 */
bool ScaffoldGraphIntegrity( const digraphE<CLinkBundle> &graph,
			     ostream *log = 0 );

#endif

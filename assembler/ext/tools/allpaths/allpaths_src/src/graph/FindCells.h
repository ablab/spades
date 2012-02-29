/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef FIND_CELLS_H
#define FIND_CELLS_H

#include "CoreTools.h"
#include "graph/Digraph.h"

// For purposes of this code, define a cell in a directed graph to be a list of 
// vertices v1,..., vn (where the order of v2,...,vn-1 does not matter) such that:
// * edges to or from v2,...,vn-1 stay within the list;
// * edges from v1 stay within the list;
// * edges to vn stay within the list;
// * every vi is between v1 and vn;
// * there is no cycle in v1,...,vn;
//   (included in this: loops v1 --> v1 and vn --> vn are disallowed;
//    we might want to change this condition);
// * n >= 2;
// * if n = 2 there are at least two edges between v1 and v2.
// A cell is minimal if it does not contain another cell as a proper 
// subset.  The size of a cell is n.

void FindCells( const digraph& G, const size_t max_cell_size, vec< vec<int> >& cells );

#endif

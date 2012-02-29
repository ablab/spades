///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/Sepdev.h"

template digraphE<sepdev>::digraphE(digraphE<sepdev> const&, vec<int> const&);
template void digraphE<sepdev>::ToLeft ( vec<int>& to_left ) const;
template void digraphE<sepdev>::ToRight( vec<int>& to_right ) const;
template void digraphE<sepdev>::ComponentEdges( vec< vec<int> >& edges ) const;
template void digraphE<sepdev>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to, const vec<sepdev>& edges, const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj );

template digraphE<fsepdev>::digraphE(digraphE<fsepdev> const&, vec<int> const&);
template void digraphE<fsepdev>::AddVertices( int n_vertices );
template void digraphE<fsepdev>::AddEdge( const int v, const int w, const fsepdev & E );
template void digraphE<fsepdev>::ToLeft ( vec<int>& to_left ) const;
template void digraphE<fsepdev>::ToRight( vec<int>& to_right ) const;
template void digraphE<fsepdev>::ComponentEdges( vec< vec<int> >& edges ) const;

void BinaryWrite( int fd, const sepdev& sd ) {
  WriteBytes( fd, &sd, sizeof( sd ) );
}

void BinaryRead( int fd, sepdev& sd ) {
  ReadBytes( fd, &sd, sizeof( sd ) );
}

void BinaryWrite( int fd, const fsepdev& fsd ) {
  WriteBytes( fd, &fsd, sizeof( fsd ) );
}

void BinaryRead( int fd, fsepdev& fsd ) {
  ReadBytes( fd, &fsd, sizeof( fsd ) );
}

template void BinaryWrite( int, const digraphE<sepdev>& );
template void BinaryWrite( int, const digraphE<fsepdev>& );

template void BinaryRead( int, digraphE<sepdev>& );
template void BinaryRead( int, digraphE<fsepdev>& );

template const Tsepdev<double>& digraphE<Tsepdev<double> >::EdgeObject(int) const;


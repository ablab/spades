#ifndef _graph_h_
#define _graph_h_

#include <vector>
#include <ostream>
#include <cassert>
#include "../sequence.hpp"
#include "parameters.hpp"

using namespace std;

class CEdge;

class CVertex
{
	Sequence * kmer;
	typedef vector < CEdge > CEdgeArray;
	CEdgeArray edgeArray;
public:
	CVertex () {};
	CVertex ( Sequence * kmer ) : kmer (kmer) {};
	CVertex (const CVertex &cv ) {}; // TODO
	void AddEdge ( CEdge e );
};

class CEdge
{
//	Seq<MPSIZE> r; TODO
//	int pos1, pos2;
	CVertex* endVertex;
public:
	CEdge ( CVertex * endVertex, Seq<MPSIZE> * r, int pos1, int pos2 ) :
		endVertex (endVertex) {};
//	Sequence label;
};

/// A straightforward graph implementation.
class CGraph
{
	typedef vector < CVertex > CVertexArray;
	CVertexArray vertexArray;
public:
	CGraph () {}
	void AddVertex ( CVertex v );
};

extern ostream & operator << ( ostream & os, const CGraph & g );

#endif // _graph_h_

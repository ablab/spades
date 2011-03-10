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
	int hits_;
	CVertex () {};
	CVertex ( Sequence * kmer ) : kmer (kmer) {};
//	CVertex (const CVertex &cv ) {};
	void AddEdge ( CEdge e );
//	int hits() {return hits_;}
};

class CEdge
{
	Seq<MPSIZE> r;
	int pos1, pos2;
	CVertex* endVertex;
public:
	CEdge ( CVertex * endVertex, Seq<MPSIZE> * r, int pos1, int pos2 ) :
		endVertex (endVertex), pos1 (pos1), pos2 (pos2) {};
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

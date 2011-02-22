#ifndef _graph_h_
#define _graph_h_

#include <vector>
#include <ostream>
#include <cassert>
#include "../seq.hpp"
#include "parameters.h"

using namespace std;

class CEdge;

class CVertex
{
	Kmer* kmer;
	typedef vector < CEdge > CEdgeArray;
public:
	CVertex ( Kmer * kmer ) : kmer (kmer) {};
};

class CEdge
{
	CVertex* endVertex;
	Sequence label;
};

/// A straightforward graph implementation.
class CGraph
{
	typedef vector < CVertex > CVertexArray;
public:
	CGraph ( unsigned number_of_vertices );
	void AddVertex ( CVertex v );
};

extern ostream & operator << ( ostream & os, const CGraph & g );

#endif // _graph_h_

#ifndef _graph_h_
#define _graph_h_

#include <vector>
#include <ostream>
#include <cassert>
#include "../seq.hpp"

using namespace std;

class CEdge;

class CVertex
{
  
  
public:
  CVertex ( unsigned id ) : id ( id ) {};
};

class CEdge
{
  typedef vector < CVertex > CVertexArray;

  CVertex * vertex;
  CVertexArray IncomingEdges;
};

/// A straightforward graph implementation.
class CGraph
{
//protected:
  typedef vector < CVertex > CVertexArray;
  //CLiterals Literals;

public:
  CGraph ( unsigned number_of_vertices );
};

extern ostream & operator << ( ostream & os, const CGraph & g );

#endif // _graph_h_


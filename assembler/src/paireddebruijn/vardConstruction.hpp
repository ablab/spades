#ifndef VARDCONSTRUCTION_H_
#define VARDCONSTRUCTION_H_
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;

namespace vard{


Sequence* SubSeq(Sequence Seq, int direction, int CutLen = 1);
int appendLowerPath(string &edge, string &toAppend);

void clearUseOfEdgePrototypes(edgesMap &edges);
void createVertices(edgesMap &edges, PairedGraph &graph);
void createEdges(edgesMap &edges, PairedGraph &graph, bool buildEdges = true);
//int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage);
int expandDirected(edgesMap &edges, protoEdgeType &curEdge, verticesMap &verts, ll &curKmer, Sequence* &curSeq, int &EdgeCoverage, int direction);
pair<char, EdgePrototype*> findUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq , int direction, bool replace = false);
int goUniqueWay(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, pair<char, EdgePrototype*>, int &EdgeCoverage, int direction);
}


#endif /*VARDCONSTRUCTION_H_*/

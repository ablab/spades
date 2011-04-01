#ifndef VARDCONSTRUCTION_H_
#define VARDCONSTRUCTION_H_
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;
namespace vard{

void createVertices(edgesMap &edges, PairedGraph &graph);
//int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage);
int expandDirected(edgesMap &edges, verticesMap &verts, ll &startKmer, Sequence* &startSeq, int &EdgeCoverage, int direction);
int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq, int direction);
int goUniqueWay(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage, int direction);
}


#endif /*VARDCONSTRUCTION_H_*/

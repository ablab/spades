#ifndef CONSTDCONSTRUCTION_H_
#define CONSTDCONSTRUCTION_H_
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "pairedGraph.hpp"
#include "graphio.hpp"

using namespace paired_assembler;
namespace constd{

void createVertices(edgesMap &edges, PairedGraph &graph);
int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage);
int expandLeft(edgesMap &edges, verticesMap &verts, ll &startKmer, Sequence* &startSeq, int &EdgeCoverage);
int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage);
int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage);
}


#endif /*CONSTDCONSTRUCTION_H_*/

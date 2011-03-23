#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph(PairedGraph &graph);
edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired = true);
void createVertices(gvis::GraphPrinter<int> &g, edgesMap &edges,
		verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph);
void createVertices(edgesMap &edges, PairedGraph &graph);
int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandLeft(edgesMap &edges, verticesMap &verts, ll &startKmer, Sequence* &startSeq);
int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
//void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph);
//void outputLongEdges(longEdgesMap &longEdges);
//void traceReads(verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph);



#endif /*GRAPHCONSTRUCTION_H_*/

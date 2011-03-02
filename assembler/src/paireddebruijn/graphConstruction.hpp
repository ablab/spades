#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"
#include "../graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph();
edgesMap sequencesToMap(string parsed_k_sequence);
void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges);
int expandDown(edgesMap &edges, vertecesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandUp(edgesMap &edges, vertecesMap &verts, ll &startKmer, Sequence* &startSeq);
int CheckUnuqueWayUp(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayUp(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int CheckUnuqueWayDown(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayDown(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int storeVertex(vertecesMap &verts, ll newKmer, Sequence* newSeq);




#endif /*GRAPHCONSTRUCTION_H_*/

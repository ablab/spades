#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"

using namespace paired_assembler;

void constructGraph();
edgesMap sequencesToMap(string parsed_k_sequence);
void createVertices(Graph *g, edgesMap &edges);
void expandDown(edgesMap &edges, vertecesMap &verts, ll finishKmer, Sequence *finishSeq);
int CheckUnuqueWayUp(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayUp(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int CheckUnuqueWayDown(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayDown(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);





#endif /*GRAPHCONSTRUCTION_H_*/

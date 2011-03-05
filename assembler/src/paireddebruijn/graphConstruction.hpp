#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph();
edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired = true);
void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges);
int expandRight(edgesMap &edges, vertecesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandLeft(edgesMap &edges, vertecesMap &verts, ll &startKmer, Sequence* &startSeq);
int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int storeVertex(gvis::GraphScheme<int> &g, vertecesMap &verts, ll newKmer, Sequence* newSeq);




#endif /*GRAPHCONSTRUCTION_H_*/

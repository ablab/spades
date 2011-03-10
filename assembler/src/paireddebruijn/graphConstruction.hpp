#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"
#include "../graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph();
<<<<<<< HEAD
edgesMap sequencesToMap(string parsed_k_sequence);
void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges);
int expandDown(edgesMap &edges, vertecesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandUp(edgesMap &edges, vertecesMap &verts, ll &startKmer, Sequence* &startSeq);
int CheckUnuqueWayUp(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayUp(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int CheckUnuqueWayDown(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int GoUnuqueWayDown(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int storeVertex(vertecesMap &verts, ll newKmer, Sequence* newSeq);


=======
edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired = true);
void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges, verticesMap &verts, longEdgesMap &longEdges);
int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandLeft(edgesMap &edges, verticesMap &verts, ll &startKmer, Sequence* &startSeq);
int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int storeVertex(gvis::GraphScheme<int> &g, verticesMap &verts, ll newKmer, Sequence* newSeq);
void resetVertexCount();
void expandDefinite(verticesMap &verts, longEdgesMap &longEdges);
void outputLongEdges(longEdgesMap &longEdges);
void traceReads(verticesMap &verts, longEdgesMap &longEdges);
>>>>>>> c9e9a4c50cceaea6b9b70c96835b90cb2ee4055d


#endif /*GRAPHCONSTRUCTION_H_*/

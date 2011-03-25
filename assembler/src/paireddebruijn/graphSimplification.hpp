#ifndef GRAPHSIMPLIFICATION_H_
#define GRAPHSIMPLIFICATION_H_



#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount, bool NotExpandBeyondDefinite = false);
void extractDefinite(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount, int dir);
bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount);
pair<int, int> vertexDist(longEdgesMap &longEdges, PairedGraph &Graph, int vertexId);
bool isPath(Edge &e1, Edge &e2);

#endif /* GRAPHSIMPLIFICATION_H_ */

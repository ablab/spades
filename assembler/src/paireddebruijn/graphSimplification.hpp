#ifndef GRAPHSIMPLIFICATION_H_
#define GRAPHSIMPLIFICATION_H_



#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount);
bool isPath(Edge &e1, Edge &e2);

#endif /* GRAPHSIMPLIFICATION_H_ */

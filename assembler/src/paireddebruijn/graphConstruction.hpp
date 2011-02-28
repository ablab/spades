#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"

using namespace paired_assembler;

void constructGraph();
edgesMap sequencesToMap(string parsed_k_sequence);
void createVertices(Graph *g, edgesMap &edges);

#endif /*GRAPHCONSTRUCTION_H_*/

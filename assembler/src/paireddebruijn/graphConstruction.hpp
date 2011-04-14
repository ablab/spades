#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph(PairedGraph &graph);
edgesMap sequencesToMap(string parsed_k_sequence);
void appendLmers(string parsed_l_mers, edgesMap &edges);

#endif /*GRAPHCONSTRUCTION_H_*/

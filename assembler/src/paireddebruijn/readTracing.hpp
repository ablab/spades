#ifndef READTRACING_HPP_
#define READTRACING_HPP_

#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

void traceReads(verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount, int &EdgeId);
bool isPath(Edge &e1, Edge &e2);
#endif /* READTRACING_HPP_ */

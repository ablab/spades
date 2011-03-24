#ifndef GRAPHSIMPLIFICATION_H_
#define GRAPHSIMPLIFICATION_H_

#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, bool NotExpandBeyondDefinite = false);
void extractDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount);
bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount);
pair<int, int> vertexDist(longEdgesMap &longEdges, PairedGraph &Graph,
		int vertexId);
bool isPath(Edge &e1, Edge &e2);

class PairThreader {
private:
	int minIntersection_;
	PairedGraph &g_;
public:
	PairThreader(PairedGraph &g, int minIntersection = 1) :
		g_(g), minIntersection_(minIntersection) {
	}
private:
	void threadLower(vector<pair<int, Edge *> > &result, Edge *currentEdge,
			int shift, Edge *start);
public:
	vector<pair<int, Edge *> > threadLower(Edge *start);
};

#endif /* GRAPHSIMPLIFICATION_H_ */

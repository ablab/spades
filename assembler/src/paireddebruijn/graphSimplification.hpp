#ifndef GRAPHSIMPLIFICATION_H_
#define GRAPHSIMPLIFICATION_H_

#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, bool NotExpandBeyondDefinite = false);
void expandObvious(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, bool NotExpandBeyondDefinite = false);
void extractDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, int dir);
void extractDefinite(PairedGraph &graph, int dir);

inline int fistDifferenceIndex(PairedGraph &graph, Edge* edge1, Edge* edge2, int direction);

void cutShortTips(PairedGraph &graph, int MaxCutLength = 5);



bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount);
pair<int, int> vertexDist(longEdgesMap &longEdges, PairedGraph &Graph,
		int vertexId);
bool isPath(Edge &e1, Edge &e2);
bool isRealPossiblePath(PairedGraph &graph, int &EdgeId, int &FollowerId);

void unglueAll(PairedGraph &graph, int direction);
void unglueAll(PairedGraph &graph);

class PairThreader {
private:
	int minIntersection_;
	PairedGraph &g_;
public:
	PairThreader(PairedGraph &g, int minIntersection = 1) :
		minIntersection_(minIntersection), g_(g) {
	}
private:
	void threadLower(vector<pair<int, Edge *> > &result, Edge *currentEdge,
			int shift, Edge *start);
public:
	vector<pair<int, Edge *> > threadLower(Edge *start);
};

#endif /* GRAPHSIMPLIFICATION_H_ */

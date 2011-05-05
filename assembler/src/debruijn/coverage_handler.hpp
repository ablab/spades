#ifndef COVERAGE_HANDLER_HPP_
#define COVERAGE_HANDLER_HPP_

#include "utils.hpp"

namespace de_bruijn {
template<class Graph>
class CoverageHandler: public GraphActionHandler<Graph> {
private:
	Graph &graph_;

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	size_t KPlusOneMerCoverage(EdgeId edge) const {
		return (size_t) (graph_.coverage(edge) * graph_.length(edge));
	}

public:
	CoverageHandler(Graph &graph) :
		graph_(graph) {
	}

	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		size_t coverage = 0;
		for (typename vector<EdgeId>::iterator it = oldEdges.begin(); it != oldEdges.end(); ++it) {
			coverage += KPlusOneMerCoverage(*it);
		}
//		cout << "single merge coverage" << endl;
//		cout << graph_.EdgeNucls(newEdge) << " " << coverage << endl;
		graph_.SetCoverage(newEdge, coverage);
	}

	virtual void HandleGlue(EdgeId oldEdge, EdgeId newEdge) {
		graph_.IncCoverage(newEdge, KPlusOneMerCoverage(oldEdge));
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		size_t length1 = graph_.length(newEdge1);
		size_t length = graph_.length(oldEdge);
		size_t coverage = KPlusOneMerCoverage(oldEdge);
		size_t coverage1 = coverage * length1 / length;
		if(coverage1 == 0)
			coverage1 = 1;
		size_t coverage2 = coverage - coverage1;
		if(coverage2 == 0)
			coverage2 = 1;
		graph_.SetCoverage(newEdge1, coverage1);
		graph_.SetCoverage(newEdge2, coverage2);
	}

	virtual ~CoverageHandler() {

	}
};
}

#endif /* COVERAGE_HANDLER_HPP_ */

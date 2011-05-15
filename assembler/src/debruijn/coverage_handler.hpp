#ifndef COVERAGE_HANDLER_HPP_
#define COVERAGE_HANDLER_HPP_

#include "utils.hpp"

namespace de_bruijn {

template<class Graph>
class CoverageHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

private:
	Graph &g_;

	size_t KPlusOneMerCoverage(EdgeId edge) const {
		return (size_t) (g_.coverage(edge) * g_.length(edge));
	}

	template <size_t k>
	void processRead(const de_bruijn::SimpleSequenceMapper<k, Graph>& threader, Read read) {
		de_bruijn::Path<EdgeId> path = threader.MapSequence(
				Sequence(read.getSequenceString()));
		if (path.sequence().size() == 0)
			return;
		const vector<EdgeId> &sequence = path.sequence();
		for (typename vector<EdgeId>::const_iterator it = sequence.begin(); it
				!= path.sequence().end(); ++it) {
			g_.IncCoverage(*it, g_.length(*it));
		}
		g_.IncCoverage(sequence[0], -path.start_pos());
		EdgeId last = sequence[sequence.size() - 1];
		g_.IncCoverage(last, path.end_pos() - g_.length(last));
	}

public:
	CoverageHandler(Graph &g) :
		g_(g) {
		g_.AddActionHandler(this);
	}

	virtual ~CoverageHandler() {
		g_.RemoveActionHandler(this);
	}

	template<size_t k, typename Stream>
	void FillCoverage(Stream& stream, const de_bruijn::EdgeIndex<k + 1, Graph>& index) {
		de_bruijn::SimpleSequenceMapper<k, Graph> threader(g_, index);
		while (!stream.eof()) {
			Read read;
			stream >> read;
			processRead(threader, read);
		}

		cout << "Here3 coverage:" << g_.coverage(index.get(Seq<k+1>("CCAC")).first) << "   length:"<< g_.length(index.get(Seq<k+1>("CCAC")).first)<< endl;
	}

	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		size_t coverage = 0;
		for (typename vector<EdgeId>::iterator it = oldEdges.begin(); it
				!= oldEdges.end(); ++it) {
			coverage += KPlusOneMerCoverage(*it);
		}
		g_.SetCoverage(newEdge, coverage);
	}

	virtual void HandleGlue(EdgeId oldEdge, EdgeId newEdge) {
		g_.IncCoverage(newEdge, KPlusOneMerCoverage(oldEdge));
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		size_t length1 = g_.length(newEdge1);
		size_t length = g_.length(oldEdge);
		size_t coverage = KPlusOneMerCoverage(oldEdge);
		size_t coverage1 = coverage * length1 / length;
		if (coverage1 == 0)
			coverage1 = 1;
		size_t coverage2 = coverage - coverage1;
		if (coverage2 == 0)
			coverage2 = 1;
		g_.SetCoverage(newEdge1, coverage1);
		g_.SetCoverage(newEdge2, coverage2);
	}
};

}

#endif /* COVERAGE_HANDLER_HPP_ */

#ifndef PAIR_INFO_FILTERS_HPP_
#define PAIR_INFO_FILTERS_HPP_

namespace omnigraph {

template<class Graph>
class AbstractPairInfoFilter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
protected:
	virtual bool Check(PairInfo<EdgeId>) const {
		return true;
	}

	virtual bool Check(EdgeId, EdgeId) const {
		return true;
	}

private:
	const Graph& graph_;

public:
	AbstractPairInfoFilter(const Graph& graph) :
			graph_(graph) {
	}

	void Filter(const PairedInfoIndex<Graph>& index,
			PairedInfoIndex<Graph>& new_index) const {
		for (auto it = index.begin(); it != index.end(); ++it) {
			auto infos = *it;
			EdgeId edge1 = infos[0].first;
			EdgeId edge2 = infos[0].second;
			if (Check(edge1, edge2)) {
				for (auto pair_info_it = infos.begin();
						pair_info_it != infos.end(); ++pair_info_it) {
					if (Check(*pair_info_it)) {
						new_index.AddPairInfo(*pair_info_it, false);
					}
				}
			}
		}
	}

	virtual ~AbstractPairInfoFilter() {
	}
};

template<class Graph>
class JumpingPairInfoChecker: public AbstractPairInfoFilter<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const GraphDistanceFinder<Graph> finder_;

public:
	JumpingPairInfoChecker(const Graph& graph, size_t short_is,
			size_t read_length, size_t delta) :
			AbstractPairInfoFilter<Graph>(graph), finder_(graph,
					short_is, read_length, delta) {
	}

protected:
	virtual bool Check(EdgeId edge1, EdgeId edge2) {
		vector<size_t> result = finder_.GetGraphDistances(edge2, edge1);
		return result.size() > 0;
	}
};

template<class Graph>
class PairInfoWeightFilter: public AbstractPairInfoFilter<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	double weight_threshold_;

public:
	PairInfoWeightFilter(const Graph &graph, double weight_threshold) :
		AbstractPairInfoFilter<Graph>(graph), weight_threshold_(weight_threshold) {
	}

protected:
	virtual bool Check(PairInfo<EdgeId> info) {
		return math::ge(info.weight, weight_threshold_);
	}
};

}

#endif /* PAIR_INFO_FILTERS_HPP_ */

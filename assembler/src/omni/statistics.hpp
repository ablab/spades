#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

namespace omnigraph {

class AbstractStatCounter {
public:
	AbstractStatCounter() {
	}

	virtual ~AbstractStatCounter() {
	}

	virtual void Count() = 0;
//protected:
//	DECL_LOGGER("StatCounter")
};

class StatList: AbstractStatCounter {
private:
	vector<AbstractStatCounter*> to_count_;
public:
	StatList(
			vector<AbstractStatCounter*> to_count =
					vector<AbstractStatCounter*> ()) :
		to_count_(to_count) {
	}

	virtual ~StatList() {
	}

	void AddStat(AbstractStatCounter* new_stat) {
		to_count_.push_back(new_stat);
	}

	const vector<AbstractStatCounter*> stats() {
		return to_count_;
	}

	virtual void Count() {
		for (size_t i = 0; i < to_count_.size(); i++) {
			to_count_[i]->Count();
		}
	}
};

template<class Graph>
class VertexEdgeStat: public AbstractStatCounter {
private:
	Graph &graph_;
public:
	VertexEdgeStat(Graph &graph) :
		graph_(graph) {
	}

	virtual ~VertexEdgeStat() {
	}

	virtual void Count() {
		size_t edgeNumber = 0;
		size_t sum_edge_length = 0;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator) {
			edgeNumber++;
			if (graph_.coverage(*iterator) > 30) {
				sum_edge_length += graph_.length(*iterator);
			}
		}
		INFO("Vertex count=" << graph_.size() << "; Edge count=" << edgeNumber);
		INFO(
				"sum length of edges(coverage > 30, both strands)"
						<< sum_edge_length);
	}
};

template<class Graph>
class BlackEdgesStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	Path<EdgeId> path1_;
	Path<EdgeId> path2_;
public:
	BlackEdgesStat(Graph &graph, Path<EdgeId> path1, Path<EdgeId> path2) :
		graph_(graph), path1_(path1), path2_(path2) {
	}

	virtual ~BlackEdgesStat(){
	}

	virtual void Count() {
		size_t black_count = 0;
		size_t edge_count = 0;
		const vector<EdgeId> path_edges1 = path1_.sequence();
		const vector<EdgeId> path_edges2 = path2_.sequence();
		set<EdgeId> colored_edges;
		colored_edges.insert(path_edges1.begin(), path_edges1.end());
		colored_edges.insert(path_edges2.begin(), path_edges2.end());
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			edge_count++;
			if (colored_edges.count(*it) == 0) {
				black_count++;
			}
		}
		if (edge_count > 0) {
			INFO(
					"Error edges count: " << black_count << " which is "
							<< 100.0 * black_count / edge_count
							<< "% of all edges");
		} else {
			INFO(
					"Error edges count: " << black_count
							<< " which is 0% of all edges");
		}
	}
};

template<class Graph>
class NStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	Path<EdgeId> path_;
	size_t perc_;
public:
	NStat(Graph &graph, Path<EdgeId> path, size_t perc = 50) :
		graph_(graph), path_(path), perc_(perc) {
	}

	virtual ~NStat(){
	}

	virtual void Count() {
		vector<size_t> lengths;
		size_t sum_all = 0;
		for (size_t i = 0; i < path_.size(); i++) {
			lengths.push_back(graph_.length(path_[i]));
			sum_all += graph_.length(path_[i]);
		}
		sort(lengths.begin(), lengths.end());
		size_t sum = 0;
		int current = lengths.size();
		while (current > 0 && sum < perc_ * 0.01 * sum_all) {
			current--;
			sum += lengths[current];
		}
		INFO("N" << perc_ << ": " << lengths[current]);
	}
};

template<class Graph>
class SelfComplementStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
public:
	SelfComplementStat(Graph &graph) :
		graph_(graph) {
	}

	virtual ~SelfComplementStat() {
	}

	virtual void Count() {
		size_t sc_number = 0;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator)
			if (graph_.conjugate(*iterator) == (*iterator))
				sc_number++;
		INFO("Self-complement count="<< sc_number);
	}
};

}

#endif /* STATISTICS_HPP_ */

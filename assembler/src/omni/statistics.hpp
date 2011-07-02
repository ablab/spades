#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include "omni_tools.hpp"
#include <map>
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

	virtual ~BlackEdgesStat() {
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

	virtual ~NStat() {
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
		INFO("Self-complement count=" << sc_number);
	}
};

template<class Graph>
class EdgePairStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	PairedInfoIndex<Graph> &pair_info_;
public:
	EdgePairStat(Graph &graph, PairedInfoIndex<Graph> &pair_info) :
		graph_(graph), pair_info_(pair_info) {
	}

	virtual ~EdgePairStat() {
	}

	void GetPairInfo(map<pair<EdgeId, EdgeId> , double> &edge_pairs) {
		for (auto iterator = pair_info_.begin(); iterator != pair_info_.end(); ++iterator) {
			vector<PairInfo<EdgeId>> v = *iterator;
			size_t w = 0;
			for (size_t i = 0; i < v.size(); i++) {
				w += v[i].weight;
			}
			edge_pairs.insert(make_pair(make_pair(v[0].first, v[0].second), w));
		}
	}

	double CountBound(map<pair<EdgeId, EdgeId> , double> &edge_pairs) {
		double sum = 0;
		vector<double> weights;
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
			sum += iterator->second;
			weights.push_back(iterator->second);
		}
		sort(weights.begin(), weights.end());
		double sum_bound = sum * 0.0001;
		size_t bound_pos = 0;
		sum = 0;
		while (bound_pos + 1 < weights.size() && sum < sum_bound)
			bound_pos++;
		return weights[bound_pos];
	}

	void CountEdgePairs(map<pair<EdgeId, EdgeId> , double> edge_pairs,
			size_t bound) {
		size_t result;
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
//			cout << iterator->second << endl;
			if (iterator->second >= bound) {
				result++;
			}
		}
		INFO("Number of edge pairs connected with paired info: " << result);
	}

	void CountTrivialEdgePairs(map<pair<EdgeId, EdgeId> , double> edge_pairs,
			size_t bound) {
		size_t result = 0;
		TrivialEdgePairChecker<Graph> checker(graph_);
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
			if (iterator->second >= bound && !checker.Check(
					iterator->first.first, iterator->first.second)) {
				result++;
			}
		}
		INFO("Number of nontrivial edge pairs connected with paired info: " << result);
	}

	virtual void Count() {
		map<pair<EdgeId, EdgeId> , double> edge_pairs;
		GetPairInfo(edge_pairs);
//		size_t bound = CountBound(edge_pairs);
		size_t bound = 20;
		CountEdgePairs(edge_pairs, bound);
		CountTrivialEdgePairs(edge_pairs, bound);
		//		double sum = 0;
		//		vector<size_t> weights;
		//		for(auto iterator = pair_info_.begin(); iterator != pair_info_.end(); ++iterator) {
		//			vector<PairInfo<EdgeId>> v = *iterator;
		//			size_t w = 0;
		//			for(size_t i = 0; i < v.size(); i++) {
		//				w += v[i].weight();
		//			}
		//			sum += w;
		//			weights.push_back(w);
		//		}
		//		sort(weights.begin(), weights.end());
		//		double bound = sum / 100;
		//		sum = 0;
		//		size_t result = 0;
		//		TrivialEdgePairChecker<Graph> checker(graph_);
		//		for(size_t i = 0; i < weights.size(); i++) {
		//			sum += weights[i];
		//			if(sum >= bound)
		//				result++;
		//		}
		//		INFO("Number of edge-pairs: " << result);
	}
};

template<class Graph>
class UniquePathStat: public omnigraph::AbstractStatCounter {

	typedef typename Graph::EdgeId EdgeId;
	Graph& g_;
	PairedInfoIndex<Graph>& pair_info_;
	size_t insert_size_;
	size_t max_read_length_;
	double variance_delta_;

	//todo get rid of this parameter
	double weight_threshold_;

	size_t considered_edge_pair_cnt_;
	size_t unique_distance_cnt_;
	size_t non_unique_distance_cnt_;


	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
		double s = 0.0;
		for (auto it = infos.begin(); it!=infos.end(); ++it) {
			if ((*it).d > g_.length((*it).first)) {
				s += (*it).weight;
			}
		}
		return s > weight_threshold_;
	}
//	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
//		for (auto it = infos.begin(); it!=infos.end(); ++it) {
//			if ((*it).d() > g_.length((*it).first())) {
//				return true;
//			}
//		}
//		return false;
//	}
public:

	UniquePathStat(Graph& g, PairedInfoIndex<Graph>& pair_info, size_t insert_size, size_t max_read_length
			, double variance_delta, double weight_threshold)
	: g_(g)
	, pair_info_(pair_info)
	, insert_size_(insert_size)
	, max_read_length_(max_read_length)
	, variance_delta_(variance_delta)
	, weight_threshold_(weight_threshold)
	, considered_edge_pair_cnt_(0)
	, unique_distance_cnt_(0)
	, non_unique_distance_cnt_(0) {

	}

	virtual void Count() {
		for (auto it = pair_info_.begin(); it != pair_info_.end(); ++it) {
			if (ContainsPositiveDistance(*it)) {
				considered_edge_pair_cnt_++;
				PairInfo<EdgeId> delegate = (*it)[0];
				EdgeId e1 = delegate.first;
				EdgeId e2 = delegate.second;
				PathCounter<Graph> counter;
				int lower_bound = insert_size_ - 2 * max_read_length_ - g_.length(e1) - g_.length(e2);
//				cout << "IS " << insert_size_ << endl;
//				cout << "MRL " << max_read_length_ << endl;
//				cout << "Raw Lower bound " << lower_bound << endl;
//				cout << "Var delta " << variance_delta_ << endl;
//				cout << "Lower bound " << (1 - variance_delta_) * lower_bound << endl;
				PathProcessor<Graph> path_processor(g_, (1 - variance_delta_) * lower_bound
						, (1 + variance_delta_) * insert_size_, g_.EdgeEnd(e1), g_.EdgeStart(e2), counter);
				path_processor.Process();
				if (counter.count() == 1) {
					unique_distance_cnt_++;
				}
				if (counter.count() > 1) {
					non_unique_distance_cnt_++;
				}
			}
		}
		INFO("Considered " << considered_edge_pair_cnt_ << " edge pairs")
		INFO(unique_distance_cnt_ << " edge pairs connected with unique path of appropriate length")
		INFO(non_unique_distance_cnt_ << " edge pairs connected with non-unique path of appropriate length")
	}

	size_t considered_edge_pair_count() {
		return considered_edge_pair_cnt_;
	}

	size_t unique_distance_count() {
		return unique_distance_cnt_;
	}

	size_t non_unique_distance_count() {
		return non_unique_distance_cnt_;
	}
private:
	DECL_LOGGER("UniquePathStat")
};

}

#endif /* STATISTICS_HPP_ */

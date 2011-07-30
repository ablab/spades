#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include "omni_tools.hpp"
#include "xmath.h"
#include "paired_info.hpp"
#include <map>

namespace omnigraph {
using namespace math;

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

	void DeleteStats() {
		for(size_t i = 0; i < to_count_.size(); i++)
			delete to_count_[i];
		to_count_.clear();
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
//		INFO("Self-complement count failed!!! ");
		INFO("Self-complement count=" << sc_number);
	}
};

template<class Graph>
class EdgePairStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	PairedInfoIndex<Graph> &pair_info_;
	const string &output_folder_;
public:
	EdgePairStat(Graph &graph, PairedInfoIndex<Graph> &pair_info, const string &output_folder) :
		graph_(graph), pair_info_(pair_info), output_folder_(output_folder) {
	}

	virtual ~EdgePairStat() {
	}

private:
	vector<double> GetWeights(map<pair<EdgeId, EdgeId> , double> &edge_pairs ) {
		vector<double> weights;
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
			weights.push_back(iterator->second);
		}
		sort(weights.begin(), weights.end());
		return weights;
	}

	void GetPairInfo(map<pair<EdgeId, EdgeId> , double> &edge_pairs, PairedInfoIndex<Graph> &index) {
		for (auto iterator = index.begin(); iterator != index.end(); ++iterator) {
			vector<PairInfo<EdgeId>> v = *iterator;
			size_t w = 0;
			for (size_t i = 0; i < v.size(); i++) {
				w += v[i].weight;
			}
			edge_pairs.insert(make_pair(make_pair(v[0].first, v[0].second), w));
		}
	}

	void RemoveTrivial(map<pair<EdgeId, EdgeId> , double> &edge_pairs) {
		TrivialEdgePairChecker<Graph> checker(graph_);
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();) {
			if (checker.Check(iterator->first.first, iterator->first.second)) {
				edge_pairs.erase(iterator++);
			} else {
				++iterator;
			}
		}
	}

//	void RemoveUntrustful(map<pair<EdgeId, EdgeId> , double> &edge_pairs, double bound) {
//		vector<double> weights;
//		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
//			weights.push_back(iterator->second);
//		}
//		sort(weights.begin(), weights.end());
//
//		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();) {
//			if(iterator->second < bound) {
//				edge_pairs.erase(iterator++);
//			} else {
//				++iterator;
//			}
//		}
//	}

public:
	vector<pair<int, double>> ComulativeHistogram(vector<double> weights) {
		vector<pair<int, double>> result;
		size_t cur = weights.size() - 1;
		size_t max = 1000;
		vector<double> res(max);
		for(int i = max - 1; i >= 0; i--) {
			while(cur >= 0 && weights[cur] >= i + 1) {
				cur--;
			}
			res[i] = weights.size() - 1 - cur;
		}
		for(size_t i = 0; i < weights.size(); i++) {
			result.push_back(make_pair(i + 1, res[i]));
		}
		return result;
	}

//	void OutputWeights(vector<double> weights, string file_name) {
//		ofstream os(file_name);
//		size_t cur = weights.size() - 1;
//		size_t max = 1000;
//		vector<double> res(max);
//		for(int i = max - 1; i >= 0; i--) {
//			while(cur >= 0 && weights[cur] >= i + 1) {
//				cur--;
//			}
//			res[i] = weights.size() - 1 - cur;
//		}
//		for(size_t i = 0; i < weights.size(); i++) {
//			os << i + 1 << " " << res[i] << endl;
//		}
//		os.close();
//	}

	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
		for (auto it = infos.begin(); it!=infos.end(); ++it) {
			if ((*it).d > graph_.length((*it).first)) {
				return true;
			}
		}
		return false;
	}

	virtual void Count() {
//		OutputWeights(GetWeights(edge_pairs), output_folder_ + "pair_info_weights.txt");
		PairedInfoIndex<Graph> new_index(graph_);
		PairInfoFilter<Graph>(graph_, 40).Filter(pair_info_, new_index);
//		RemoveUntrustful(edge_pairs, 40);
		map<pair<EdgeId, EdgeId> , double> edge_pairs;
		TrivialEdgePairChecker<Graph> checker(graph_);
		size_t nontrivial = 0;
		size_t pair_number = 0;
		for(auto iterator = new_index.begin(); iterator != new_index.end(); ++iterator) {
			vector<PairInfo<EdgeId>> info = *iterator;
			if(ContainsPositiveDistance(info)) {
				pair_number++;
				if(checker.Check(info[0].first, info[0].second)) {
					nontrivial++;
				}
			}
		}
		GetPairInfo(edge_pairs, new_index);
		INFO("Number of edge pairs connected with paired info: " << pair_number);
		RemoveTrivial(edge_pairs);
		INFO("Number of nontrivial edge pairs connected with paired info: " << nontrivial);
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

	bool draw_pictures_;

//	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
//		double s = 0.0;
//		for (auto it = infos.begin(); it!=infos.end(); ++it) {
//			if ((*it).d > g_.length((*it).first)) {
//				s += (*it).weight;
//			}
//		}
//		return s > weight_threshold_;
//	}
	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
		for (auto it = infos.begin(); it!=infos.end(); ++it) {
			if ((*it).d > g_.length((*it).first)) {
				return true;
			}
		}
		return false;
	}
public:

	UniquePathStat(Graph& g, PairedInfoIndex<Graph>& pair_info, size_t insert_size, size_t max_read_length
			, double variance_delta, double weight_threshold, bool draw_pictures = false)
	: g_(g)
	, pair_info_(pair_info)
	, insert_size_(insert_size)
	, max_read_length_(max_read_length)
	, variance_delta_(variance_delta)
	, weight_threshold_(weight_threshold)
	, considered_edge_pair_cnt_(0)
	, unique_distance_cnt_(0)
	, non_unique_distance_cnt_(0)
	, draw_pictures_(draw_pictures)
	{

	}

	virtual void Count() {
		PairedInfoIndex<Graph> filtered_index(g_);
		PairInfoFilter<Graph>(g_, 40).Filter(pair_info_, filtered_index);

		for (auto it = filtered_index.begin(); it != filtered_index.end(); ++it) {
			if (ContainsPositiveDistance(*it)) {
				considered_edge_pair_cnt_++;
				PairInfo<EdgeId> delegate = (*it)[0];
				EdgeId e1 = delegate.first;
				EdgeId e2 = delegate.second;

//				cout << "Finding paths between edges " << e1 << " and " << e2 << endl;
				NonEmptyPathCounter<Graph> counter(g_);
//				VertexLablerCallback<Graph> graph_labeler(g_);
//				CompositeCallback<Graph> composite_callback;
//				composite_callback.AddProcessor(counter);
//				composite_callback.AddProcessor(graph_labeler);
				//todo delete unnecessary parentheses and casts
				int lower_bound = ((int) insert_size_) - 2 * ((int) max_read_length_) - ((int) g_.length(e1)) - ((int) g_.length(e2));
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

template<class Graph>
class EstimationQualityStat: public omnigraph::AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef PairInfo<EdgeId> Info;
	typedef typename PairedInfoIndex<Graph>::PairInfos Infos;
	PairedInfoIndex<Graph>& estimated_pair_info_;
	PairedInfoIndex<Graph>& etalon_pair_info_;
	vector<double> false_positive_weights_;
	vector<double> perfect_match_weights_;
	//todo make count fair (now all the cluster weights are counted several times)
	vector<double> imperfect_match_weights_;
	size_t false_negative_count_;
private:

	void HandleFalsePositive(const Info& estimated) {
		false_positive_weights_.push_back(estimated.weight);
	}

	void HandleFalseNegative(const Info& etalon) {
		false_negative_count_++;
	}

	void HandlePairsNotInEtalon(const set<pair<EdgeId, EdgeId>>& pairs_in_etalon) {
		for (auto it = estimated_pair_info_.begin(); it != estimated_pair_info_.end(); ++it) {
			Infos estimated_infos = *it;
			EdgeId first = estimated_infos[0].first;
			EdgeId second = estimated_infos[0].second;
			if (pairs_in_etalon.count(make_pair(first, second)) == 0) {
				for_each(estimated_infos.begin(), estimated_infos.end(), HandleFalsePositive);
			}
		}
	}

	bool InfoLess(const Info& a, const Info& b) {
		return le(a.d + a.variance, b.d - b.variance);
	}

	bool PerfectMatch(const Info& etalon, const Info& estimated) {
		if (eq(etalon.d, estimated.d) && eq(estimated.variance, 0.)) {
			perfect_match_weights_.push_back(estimated.weight);
			return true;
		}
		return false;
	}

	bool ImperfectMatch(const Info& etalon, const Info& estimated) {
		if (ge(etalon.d > estimated.d - estimated.variance) && le(etalon.d < estimated.d + estimated.variance)) {
			imperfect_match_weights_.push_back(estimated.weight);
			return true;
		}
		return false;
	}

	void ProcessInfos(const Infos& etalon_infos, const Infos& estimated_infos) {
		size_t estimated_idx = 0;
		bool last_matched = false;
		for (size_t etalon_idx = 0; etalon_idx < etalon_infos.size(); ++etalon_idx) {
			Info etalon_info = etalon_infos[etalon_idx];
			while (estimated_idx < estimated_infos.size() && InfoLess(estimated_infos[estimated_idx], etalon_info)) {
				HandleFalsePositive(estimated_infos[estimated_idx]);
				estimated_idx++;
			}
			if (estimated_idx != estimated_infos.size()
					&& (PerfectMatch(etalon_info, estimated_infos[estimated_idx])
							|| ImperfectMatch(etalon_info, estimated_infos[estimated_idx]))) {
				last_matched = true;
			} else {
				HandleFalseNegative(etalon_info);
			}
		}
		if (last_matched)
			estimated_idx++;
		while (estimated_idx < estimated_infos.size()) {
			HandleFalsePositive(estimated_infos[estimated_idx]);
		}
	}

public:
	EstimationQualityStat(PairedInfoIndex<Graph>& estimated_pair_info, PairedInfoIndex<Graph>& etalon_pair_info) :
		estimated_pair_info_(estimated_pair_info), etalon_pair_info_(etalon_pair_info), false_negative_count_(0) {
	}

	virtual ~EstimationQualityStat() {
	}

	virtual void Count() {
		set<pair<EdgeId, EdgeId>> pairs_in_etalon;
		for (auto it = etalon_pair_info_.begin(); it != etalon_pair_info_.end(); ++it) {
			Infos etalon_infos = *it;
			EdgeId first = etalon_infos[0].first;
			EdgeId second = etalon_infos[0].second;
			pairs_in_etalon.insert(first, second);

			Infos estimated_infos = estimated_pair_info_.GetEdgePairInfo(first, second);
			ProcessInfos(etalon_infos, estimated_infos);
		}
		HandlePairsNotInEtalon(pairs_in_etalon);
	}

	vector<double> false_positive_weights() {
		sort(false_positive_weights_.begin(), false_positive_weights_.end());
		return false_positive_weights_;
	}
	vector<double> perfect_match_weights() {
		sort(perfect_match_weights_.begin(), perfect_match_weights_.end());
		return perfect_match_weights_;
	}
	//todo make count fair (now all the cluster weights are counted several times)
	vector<double> imperfect_match_weights() {
		sort(imperfect_match_weights_.begin(), imperfect_match_weights_.end());
		return imperfect_match_weights_;
	}
	size_t false_negative_count() {
		return false_negative_count_;
	}

};

}

#endif /* STATISTICS_HPP_ */

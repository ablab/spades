//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include "standard.hpp"
#include "omni/omni_tools.hpp"
#include "simple_tools.hpp"
#include "xmath.h"
#include "omni/paired_info.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/pair_info_filters.hpp"

#include "graphio.hpp"
#include <iostream>
#include <fstream>
#include <map>

#include "path_set.hpp" 
#include "matepair_transformer.hpp"
#include "path_set_stats.hpp"


namespace debruijn_graph {
using namespace math;
using namespace omnigraph;

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
	StatList(vector<AbstractStatCounter*> to_count =
			vector<AbstractStatCounter*>()) :
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
		for (size_t i = 0; i < to_count_.size(); i++)
			delete to_count_[i];
		to_count_.clear();
	}
};

template<class Graph>
class VertexEdgeStat: public AbstractStatCounter {
private:
	const Graph &graph_;
public:
	VertexEdgeStat(const Graph &graph) :
			graph_(graph) {
	}

	virtual ~VertexEdgeStat() {
	}

	size_t vertices() {
		return graph_.size();
	}

	size_t edges() {
		size_t edgeNumber = 0;
		size_t sum_edge_length = 0;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();
				++iterator) {
			edgeNumber++;
//			if (graph_.coverage(*iterator) > 30) {
				sum_edge_length += graph_.length(*iterator);
//			}
		}
		return edgeNumber;
	}

	size_t edge_length() {
		size_t sum_edge_length = 0;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();
				++iterator) {
			if (graph_.coverage(*iterator) > 30) {
				sum_edge_length += graph_.length(*iterator);
			}
		}
		return sum_edge_length;
	}

	virtual void Count() {
		INFO(
				"Vertex count=" << vertices() << "; Edge count=" << edges());
		INFO(
				"sum length of edges " << edge_length());
	}
};

template<class Graph>
class BlackEdgesStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	Path<EdgeId> path1_;
	Path<EdgeId> path2_;
public:
	BlackEdgesStat(const Graph &graph, Path<EdgeId> path1, Path<EdgeId> path2) :
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
		size_t sum_length = 0;
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			edge_count++;
			if (colored_edges.count(*it) == 0) {
				black_count++;
				sum_length += graph_.length(*it);
			}
		}
		if (edge_count > 0) {
			INFO(
					"Error edges count: " << black_count << " which is " << 100.0 * black_count / edge_count << "% of all edges");
			INFO(
					"Total length of all black edges: " << sum_length << ". While double genome length is " << (2 * cfg::get().ds.reference_genome.size()));
		} else {
			INFO(
					"Error edges count: " << black_count << " which is 0% of all edges");
		}
	}
};

template<class Graph>
class NStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	Path<EdgeId> path_;
	size_t perc_;
public:
	NStat(const Graph &graph, Path<EdgeId> path, size_t perc = 50) :
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
		size_t current = lengths.size();
		while (current > 0 && sum < perc_ * 0.01 * sum_all) {
			current--;
			sum += lengths[current];
		}
		if(current < lengths.size())
			INFO("N" << perc_ << ": " << lengths[current]);
	}
};

template<class Graph>
class IsolatedEdgesStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	set<EdgeId> black_edges_;
	vector<size_t> lengths;
public:
	IsolatedEdgesStat(const Graph &graph, Path<EdgeId> path1,
			Path<EdgeId> path2) :
			graph_(graph) {
		for (auto it = graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			black_edges_.insert(*it);
		}
		for (size_t i = 0; i < path1.size(); i++) {
			black_edges_.erase(path1[i]);
		}
		for (size_t i = 0; i < path2.size(); i++) {
			black_edges_.erase(path2[i]);
		}
	}

	virtual ~IsolatedEdgesStat() {
	}

	virtual void Count() {
		lengths.clear();
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			EdgeId edge = *it;
			if (graph_.IsDeadEnd(graph_.EdgeEnd(edge))
					&& graph_.IsDeadStart(graph_.EdgeStart(edge))
					&& black_edges_.count(edge) == 0) {
				lengths.push_back(graph_.length(edge));
			}
		}INFO("Isolated not black edges: " << lengths.size());
		WriteLengths(cfg::get().output_dir, "isolated_edges.txt");
	}

	void WriteLengths(string folder_name, string file_name) {
		ofstream os;
		os.open((folder_name + "/" + file_name).c_str());
		WriteLengths(os);
		os.close();
	}

	void WriteLengths(ostream &os) {
		sort(lengths.begin(), lengths.end());
		for (size_t i = 0; i < lengths.size(); i++) {
			os << lengths[i] << endl;
		}
	}
};

template<class Graph>
class SelfComplementStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
public:
	SelfComplementStat(const Graph &graph) :
			graph_(graph) {
	}

	virtual ~SelfComplementStat() {
	}

	virtual void Count() {
		size_t sc_number = 0;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();
				++iterator)
			if (graph_.conjugate(*iterator) == (*iterator))
				sc_number++;
		//		INFO("Self-complement count failed!!! ");
		INFO("Self-complement count=" << sc_number);
	}
};

template<class Graph>
class EdgePairStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const PairedInfoIndex<Graph> &pair_info_;
	const string &output_folder_;
public:
	EdgePairStat(const Graph &graph, const PairedInfoIndex<Graph> &pair_info,
			const string &output_folder) :
			graph_(graph), pair_info_(pair_info), output_folder_(output_folder) {
	}

	virtual ~EdgePairStat() {
	}

private:
	vector<double> GetWeights(map<pair<EdgeId, EdgeId> , double> &edge_pairs) {
		vector<double> weights;
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();
				++iterator) {
			weights.push_back(iterator->second);
		}
		sort(weights.begin(), weights.end());
		return weights;
	}

	void GetPairInfo(map<pair<EdgeId, EdgeId> , double> &edge_pairs,
	PairedInfoIndex<Graph> &index) {
		for (auto iterator = index.begin(); iterator != index.end();
				++iterator) {
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
		for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();
				) {
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
		int cur = weights.size() - 1;
		size_t max = 1000;
		vector<double> res(max);
		for (int i = max - 1; i >= 0; i--) {
			while (cur >= 0 && weights[cur] >= i + 1) {
				cur--;
			}
			res[i] = weights.size() - 1 - cur;
		}
		for (size_t i = 0; i < weights.size(); i++) {
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
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if ((*it).d > graph_.length((*it).first)) {
				return true;
			}
		}
		return false;
	}

	virtual void Count() {
		//		OutputWeights(GetWeights(edge_pairs), output_folder_ + "pair_info_weights.txt");
		PairedInfoIndex<Graph> new_index(graph_);
		PairInfoWeightFilter<Graph>(graph_, 40).Filter(pair_info_, new_index);
		//		RemoveUntrustful(edge_pairs, 40);
		map<pair<EdgeId, EdgeId> , double> edge_pairs;
		TrivialEdgePairChecker<Graph> checker(graph_);
		size_t nontrivial = 0;
		size_t pair_number = 0;
		for (auto iterator = new_index.begin(); iterator != new_index.end();
				++iterator) {
			vector<PairInfo<EdgeId>> info = *iterator;
			if (ContainsPositiveDistance(info)) {
				pair_number++;
				if (checker.Check(info[0].first, info[0].second)) {
					nontrivial++;
				}
			}
		}
		GetPairInfo(edge_pairs, new_index);
		INFO(
				"Number of edge pairs connected with paired info: " << pair_number);
		RemoveTrivial(edge_pairs);
		INFO(
				"Number of nontrivial edge pairs connected with paired info: " << nontrivial);
	}
};

template<class Graph>
class UniquePathStat: public AbstractStatCounter {

	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const PairedInfoIndex<Graph>& filtered_index_;
	size_t insert_size_;
	size_t max_read_length_;
	size_t gap_;
	double variance_delta_;

	size_t considered_edge_pair_cnt_;
	size_t unique_distance_cnt_;
	size_t non_unique_distance_cnt_;

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
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if ((*it).d > g_.length((*it).first)) {
				return true;
			}
		}
		return false;
	}
public:

	UniquePathStat(const Graph& g, const PairedInfoIndex<Graph>& filtered_index,
			size_t insert_size, size_t max_read_length, double variance_delta) :
			g_(g), filtered_index_(filtered_index), insert_size_(insert_size), max_read_length_(
					max_read_length), gap_(insert_size_ - 2 * max_read_length_), variance_delta_(
					variance_delta), considered_edge_pair_cnt_(0), unique_distance_cnt_(
					0), non_unique_distance_cnt_(0) {

	}

	virtual void Count() {
		//		PairedInfoIndex<Graph> filtered_index(g_);
		//		PairInfoFilter < Graph > (g_, 40).Filter(pair_info_, filtered_index);

		for (auto it = filtered_index_.begin(); it != filtered_index_.end();
				++it) {
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
				PathProcessor<Graph> path_processor(
						g_,
						omnigraph::PairInfoPathLengthLowerBound(g_.k(),
								g_.length(e1), g_.length(e2), gap_,
								variance_delta_),
						omnigraph::PairInfoPathLengthUpperBound(g_.k(),
								insert_size_, variance_delta_), g_.EdgeEnd(e1),
						g_.EdgeStart(e2), counter);
				path_processor.Process();
				if (counter.count() == 1) {
					unique_distance_cnt_++;
				}
				if (counter.count() > 1) {
					non_unique_distance_cnt_++;

				}
			}
		}INFO("Considered " << considered_edge_pair_cnt_ << " edge pairs")INFO(
				unique_distance_cnt_ << " edge pairs connected with unique path of appropriate length")
		INFO(
				non_unique_distance_cnt_ << " edge pairs connected with non-unique path of appropriate length")
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
class MatePairTransformStat: public AbstractStatCounter {

	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const PairedInfoIndex<Graph>& pair_info_;

	size_t considered_dist_cnt_;
	size_t unique_distance_cnt_;
	size_t non_unique_distance_cnt_;

	void ProcessInfos(const vector<PairInfo<EdgeId>>& infos) {
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			PairInfo<EdgeId> info = *it;
			if (gr(info.d, 0.)) {
				if (info.variance == 0) {
					EdgeId e1 = info.first;
					EdgeId e2 = info.second;

					PathStorageCallback<Graph> counter(g_);

					PathProcessor<Graph> path_processor(g_,
							info.d - g_.length(e1), info.d - g_.length(e1),
							g_.EdgeEnd(e1), g_.EdgeStart(e2), counter);
					path_processor.Process();
					TRACE(
							"Edges"<< e1 <<" : " << e2 << ": " << info.weight << " : " << info.d);
					TRACE("Path Numbs" << counter.size());

					if (counter.size() == 1) {
						unique_distance_cnt_++;
					}
					if (counter.size() > 1) {
						non_unique_distance_cnt_++;
					}
				} else {
					non_unique_distance_cnt_++;
				}
				considered_dist_cnt_++;
			}
		}
	}

public:

	MatePairTransformStat(const Graph& g,
			const PairedInfoIndex<Graph>& pair_info) :
			g_(g), pair_info_(pair_info), considered_dist_cnt_(0), unique_distance_cnt_(
					0), non_unique_distance_cnt_(0) {

	}

	virtual void Count() {

        /*
        PathSetIndexData<EdgeId> PII ;
        PathSetIndexData<EdgeId> PIIFilter ;
        MatePairTransformer<Graph> transformer(g_, pair_info_);
        transformer.Transform(PII);
        PathSetIndex<EdgeId> PI(PII);
        INFO("STATISTIC");
        PI.Process(PIIFilter);
        INFO("END STATISTIC");
        for(auto iter = PII.begin(); iter != PII.end() ; ++iter)
        {
            INFO( *iter);
        } 
        
//        PathSetStats<Graph>  pathsetStatistic(g_,PII , PIIFilter);
//        pathsetStatistic.Count();


        INFO("FILTERED");
        for(auto iter = PIIFilter.begin() ; iter != PIIFilter.end(); ++iter)
        {
            INFO(*iter);
        }
        INFO("Extension");
        for(auto iter = PIIFilter.begin(); iter != PIIFilter.end() ; ++iter)
        {
            PathSet<EdgeId> first = *iter;
            INFO(first);
            vector<PathSet<EdgeId>> extends;
            size_t offSet = 0;
            if(first.paths.size() > 1)
                offSet = 1;
            PI.FindExtension(PIIFilter,first, extends, offSet);
            for(size_t i = 0 ; i < extends.size() ; ++i)
            {
                if(offSet == 0)
                {
                    INFO("CAN BE EXTENDED BY " );
                }
                else
                {
                    INFO("CAN BE EXTENDED BY LONGER");
                }
                INFO(extends[i]);
            }
        }

        */

        
		for (auto it = pair_info_.begin(); it != pair_info_.end(); ++it) {
			vector<PairInfo<EdgeId>> infos = *it;
			ProcessInfos(infos);
		}INFO(
				"Considered " << considered_dist_cnt_ << " edge pair distances (including trivial)");
		INFO(
				unique_distance_cnt_ << " edge distances connected with unique path of appropriate length");
		INFO(
				non_unique_distance_cnt_ << " edge distances connected with non-unique path of appropriate length");
	}

	size_t considered_edge_pair_count() {
		return considered_dist_cnt_;
	}

	size_t unique_distance_count() {
		return unique_distance_cnt_;
	}

	size_t non_unique_distance_count() {
		return non_unique_distance_cnt_;
	}
private:
	DECL_LOGGER("MatePairTransformStat")
};

template<class Graph>
class UniqueDistanceStat: public AbstractStatCounter {
	typedef omnigraph::PairedInfoIndex<Graph> PairedIndex;

	const PairedIndex& paired_info_;
	size_t unique_;
	size_t non_unique_;
public:

	UniqueDistanceStat(const PairedIndex& paired_info) :
			paired_info_(paired_info), unique_(0), non_unique_(0) {

	}

	virtual ~UniqueDistanceStat() {

	}

	virtual void Count() {
		for (auto it = paired_info_.begin(); it != paired_info_.end(); ++it) {
			VERIFY((*it).size() > 0);
			if ((*it).size() > 1) {
				non_unique_++;
				//				for (auto info_it = (*it).begin(); info_it != (*it).end(); ++info_it) {
				//					//todo
				//				}
			} else {
				unique_++;
			}
		}INFO(unique_ << " unique edge distances");
		INFO(non_unique_ << " non unique edge distances");
	}

	size_t unique() {
		return unique_;
	}

	size_t non_unique() {
		return non_unique_;
	}
};

template<class Graph>
class EstimationQualityStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef PairInfo<EdgeId> Info;
	typedef typename PairedInfoIndex<Graph>::PairInfos Infos;
	//input fields
	const Graph &graph_;
	const IdTrackHandler<Graph>& int_ids_;
	const EdgeQuality<Graph>& quality_;
	const PairedInfoIndex<Graph>& pair_info_;
	const PairedInfoIndex<Graph>& estimated_pair_info_;
	const PairedInfoIndex<Graph>& etalon_pair_info_;

	//output fields
	PairedInfoIndex<Graph> false_positives_;
	PairedInfoIndex<Graph> perfect_matches_;
	PairedInfoIndex<Graph> imperfect_matches_;
	PairedInfoIndex<Graph> false_negatives_;

//	PairedInfoIndex<Graph> false_positive_weights_;
//	set<Info> false_positive_infos_;
//	vector<double> perfect_match_weights_;
////(weight, estimated_variance - actual_variance, number of etalon points)
//	vector<pair<pair<double, double> , size_t>> imperfect_match_stat_;
//	size_t false_negative_count_;
//	vector<Info> false_negative_infos_;

	bool CheckInterestInInfo(const Info& info) {
        if (math::ls(info.d, 0.)) return false;
        if (info.first == info.second && math::eq(info.d, 0.)) return false;
		return quality_.IsPositiveQuality(info.first)
						&& quality_.IsPositiveQuality(info.second) && math::gr(info.weight, 0.);
	}

	void HandleFalsePositive(const Info& estimated) {
//		DEBUG("Handling false positive " << estimated);
		if (CheckInterestInInfo(estimated))
			false_positives_.AddPairInfo(estimated, false);
	}

	void HandleFalseNegative(const Info& etalon) {
		if (CheckInterestInInfo(etalon))
			false_negatives_.AddPairInfo(etalon, false);
	}

	void HandlePerfectMatch(const Info& etalon, const Info& estimated) {
		if (CheckInterestInInfo(estimated))
			perfect_matches_.AddPairInfo(estimated, false);
	}

	void HandleImperfectMatch(const Info &estimated_cluster,
			const Infos& etalon_matches) {
		if (CheckInterestInInfo(estimated_cluster))
			imperfect_matches_.AddPairInfo(estimated_cluster, false);
//		double etalon_variance = etalon_matches[etalon_matches.size() - 1].d
//				- etalon_matches[0].d;
//		imperfect_match_stat_.push_back(
//				make_pair(
//						make_pair(estimated_cluster.weight,
//								estimated_cluster.variance - etalon_variance),
//						etalon_matches.size()));
	}

//	void Flush() {
//		ProcessImperfectMatch(last_estimated_imperfect_match_,
//				last_etalon_imperfect_matches_);
//	}

	void HandlePairsNotInEtalon(
			const set<pair<EdgeId, EdgeId>>& pairs_in_etalon) {
		for (auto it = estimated_pair_info_.begin();
				it != estimated_pair_info_.end(); ++it) {
			Infos estimated_infos = *it;
			EdgeId first = estimated_infos[0].first;
			EdgeId second = estimated_infos[0].second;
			if (pairs_in_etalon.count(make_pair(first, second)) == 0) {
				//				for_each(estimated_infos.begin(), estimated_infos.end(),
				//						boost::bind(&EstimationQualityStat::HandleFalsePositive, this, _1));

				for (auto it2 = estimated_infos.begin();
						it2 != estimated_infos.end(); ++it2) {
					HandleFalsePositive(*it2);
				}
			}
		}
	}

	bool InfoLess(const Info& a, const Info& b) {
		if (eq(a.variance, 0.) && eq(b.variance, 0.)) {
			return ls(a.d, b.d);
		}
		return ls(a.d + a.variance, b.d - b.variance);
	}

	bool IsPerfectMatch(const Info& etalon, const Info& estimated) {
		return le(etalon.d, estimated.d) && ge(etalon.d, estimated.d)
				&& eq(estimated.variance, 0.);
	}

	bool IsImperfectMatch(const Info& etalon, const Info& estimated) {
		return ge(etalon.d, estimated.d - estimated.variance)
				&& le(etalon.d, estimated.d + estimated.variance);
	}

	size_t Move(size_t estimated_idx, const Infos &estimated_infos) {
		estimated_idx++;
		while (estimated_idx < estimated_infos.size()
				&& math::eq(estimated_infos[estimated_idx].weight, 0.))
			estimated_idx++;
		return estimated_idx;
		return 0;
	}

	size_t InitIdx(const Infos &pair_infos) {
		return Move(-1, pair_infos);
	}

	void ProcessInfos(const Infos& etalon_infos, const Infos& estimated_infos) {
		size_t etalon_idx = InitIdx(etalon_infos);
		for (size_t estimated_idx = InitIdx(estimated_infos);
				estimated_idx < estimated_infos.size();
				estimated_idx = Move(estimated_idx, estimated_infos)) {
			while (estimated_idx < estimated_infos.size()
					&& (etalon_idx == etalon_infos.size()
							|| InfoLess(estimated_infos[estimated_idx],
									etalon_infos[etalon_idx]))) {
				HandleFalsePositive(estimated_infos[estimated_idx]);
				estimated_idx = Move(estimated_idx, estimated_infos);
			}
			if (estimated_idx == estimated_infos.size()) {
				break;
			}
			while (etalon_idx < etalon_infos.size()
					&& InfoLess(etalon_infos[etalon_idx],
							estimated_infos[estimated_idx])) {
				HandleFalseNegative(etalon_infos[etalon_idx]);
				etalon_idx = Move(etalon_idx, etalon_infos);
			}
			if (etalon_idx == etalon_infos.size()) {
				continue;
			}
			if (IsPerfectMatch(etalon_infos[etalon_idx],
					estimated_infos[estimated_idx])) {
				while (etalon_idx < etalon_infos.size()
						&& IsPerfectMatch(etalon_infos[etalon_idx],
								estimated_infos[estimated_idx])) {
					HandlePerfectMatch(etalon_infos[etalon_idx],
							estimated_infos[estimated_idx]);
					etalon_idx = Move(etalon_idx, etalon_infos);
				}
			} else {
				vector<PairInfo<EdgeId> > cluster_hits;
				while (etalon_idx < etalon_infos.size()
						&& IsImperfectMatch(etalon_infos[etalon_idx],
								estimated_infos[estimated_idx])) {
					cluster_hits.push_back(etalon_infos[etalon_idx]);
					etalon_idx = Move(etalon_idx, etalon_infos);
				}
				if (cluster_hits.size() == 0) {
					HandleFalsePositive(estimated_infos[estimated_idx]);
				} else {
					HandleImperfectMatch(estimated_infos[estimated_idx],
							cluster_hits);
				}
			}
		}
		//		for (size_t etalon_idx = 0; etalon_idx < etalon_infos.size(); ++etalon_idx) {
		//			Info etalon_info = etalon_infos[etalon_idx];
		////			cout << "here" << endl;
		//			while (estimated_idx < estimated_infos.size() && InfoLess(estimated_infos[estimated_idx], etalon_info)) {
		//				HandleFalsePositive(estimated_infos[estimated_idx]);
		//				estimated_idx++;
		////				cout << "here1" << endl;
		//			}
		////			cout << "here2" << endl;
		//			if (estimated_idx != estimated_infos.size()
		//					&& (HandleIfPerfectMatch(etalon_info, estimated_infos[estimated_idx])
		//							|| HandleIfImperfectMatch(etalon_info, estimated_infos[estimated_idx]))) {
		//				last_matched = true;
		//			} else {
		//				HandleFalseNegative(etalon_info);
		//			}
		//		}
		//		if (last_matched)
		//			estimated_idx++;
		while (etalon_idx < etalon_infos.size()) {
			//			DEBUG("Handling false positives beyond all etalons");
			HandleFalseNegative(etalon_infos[etalon_idx]);
			etalon_idx = Move(etalon_idx, etalon_infos);
		}
		//		Flush();
	}

//	void ReportFalsePositiveWeights() {
//		sort(false_positive_weights_.begin(), false_positive_weights_.end());
//
//		INFO("False positive count: " << false_positive_weights_.size());
//	}
//
//	void ReportPerfectMatchWeights() {
//		sort(perfect_match_weights_.begin(), perfect_match_weights_.end());
//		INFO("Perfect match count: " << perfect_match_weights_.size());
//	}
//
//	void ReportImperfectMatchWeights() {
//		sort(imperfect_match_stat_.begin(), imperfect_match_stat_.end());
//		//todo do something better
//		INFO("Imperfect match count: " << imperfect_match_stat_.size());
//	}
//
//	void FalseNegativeCount() {
//		INFO("False negative count: " << false_negative_count_);
//	}

public:
	EstimationQualityStat(const Graph &graph,
			const IdTrackHandler<Graph>& int_ids,
			const EdgeQuality<Graph>& quality,
			const PairedInfoIndex<Graph>& pair_info,
			const PairedInfoIndex<Graph>& estimated_pair_info,
			const PairedInfoIndex<Graph>& etalon_pair_info) :
			graph_(graph), int_ids_(int_ids), quality_(quality), pair_info_(pair_info), estimated_pair_info_(
					estimated_pair_info), etalon_pair_info_(etalon_pair_info), false_positives_(
					graph_), perfect_matches_(graph_), imperfect_matches_(
					graph_), false_negatives_(graph_) {
	}

	virtual ~EstimationQualityStat() {
	}

	virtual void Count() {
		INFO("Counting distance estimation statistics");
		set<pair<EdgeId, EdgeId>> pairs_in_etalon;
		//		DEBUG("Handling pairs present in etalon information");
		for (auto it = etalon_pair_info_.begin(); it != etalon_pair_info_.end();
				++it) {
			Infos etalon_infos = *it;
			EdgeId first = etalon_infos[0].first;
			EdgeId second = etalon_infos[0].second;
			pairs_in_etalon.insert(make_pair(first, second));

			Infos estimated_infos = estimated_pair_info_.GetEdgePairInfo(first,
					second);
			//			DEBUG("Processing distances for pair " << first << ", " << second);
			ProcessInfos(etalon_infos, estimated_infos);
		}
		//		DEBUG("Handling pairs that are not in etalon information");
		HandlePairsNotInEtalon(pairs_in_etalon);

		INFO("FPR: " << fpr());
		INFO("FNR: " << fnr());
		INFO("Distance estimation statistics counted");
	}

	const PairedInfoIndex<Graph>& false_positives() {
		return false_positives_;
	}

	const PairedInfoIndex<Graph>& perfect_matches() {
		return perfect_matches_;
	}

	const PairedInfoIndex<Graph>& imperfect_matches() {
		return imperfect_matches_;
	}

	const PairedInfoIndex<Graph>& false_negatives() {
		return false_negatives_;
	}

	double fpr() {
		return 1. * false_positives_.size() / estimated_pair_info_.size();
	}

	double fnr() {
		return 1. * false_negatives_.size() / etalon_pair_info_.size();
	}

	void SaveStats(const string& dir_name) {
		//saving results
		INFO("Saving estimation statistic");
		make_dir(dir_name);
		typename PrinterTraits<Graph>::Printer printer(graph_, int_ids_);
		printer.savePaired(dir_name + "fp", false_positives_);
		printer.savePaired(dir_name + "pm", perfect_matches_);
		printer.savePaired(dir_name + "im", imperfect_matches_);
		printer.savePaired(dir_name + "fn", false_negatives_);
		INFO("Estimation statistics saved");
	}

//	vector<double> false_positive_weights() {
//		sort(false_positive_weights_.begin(), false_positive_weights_.end());
//		return false_positive_weights_;
//	}
//	vector<double> perfect_match_weights() {
//		sort(perfect_match_weights_.begin(), perfect_match_weights_.end());
//		return perfect_match_weights_;
//	}
//
//	vector<pair<pair<double, double> , size_t>> imperfect_match_weights() {
//		sort(imperfect_match_stat_.begin(), imperfect_match_stat_.end());
//		return imperfect_match_stat_;
//	}
//
//	size_t false_negative_count() {
//		return false_negative_count_;
//	}

//	void WriteFalseNegativeGaps(const string &file_name) {
//		ofstream stream;
//		stream.open(file_name);
//		vector<double> to_print;
//		//		for (size_t i = 0; i < false_negative_infos_.size(); i++) {
//		//			if (false_negative_infos_[i].d > 0)
//		//				to_print.push_back(
//		//						false_negative_infos_[i].d - graph_.length(
//		//								false_negative_infos_[i].first));
//		//		}
//		//		sort(to_print.begin(), to_print.end());
//		//		copy(to_print.begin(), to_print.end(),
//		//				ostream_iterator<double> (stream, "\n"));
//		for (size_t i = 0; i < false_negative_infos_.size(); i++) {
//			stream << false_negative_infos_[i] << endl;
//		}
//		stream.close();
//	}
//
//	void WriteEstmationStats(const string &output_folder) {
//		ofstream stream;
//		stream.open(output_folder + "/perfect.inf");
//		copy(perfect_match_weights_.begin(), perfect_match_weights_.end(),
//				ostream_iterator<double>(stream, "\n"));
//		stream.close();
//
//		stream.open(output_folder + "/false_positive.inf");
//		copy(false_positive_weights_.begin(), false_positive_weights_.end(),
//				ostream_iterator<double>(stream, "\n"));
//		stream.close();
//		WriteWorstEdgesStat(output_folder, 1000000);
//	}

//	void WriteEdgePairInfo(const string &file_name, Infos infos) {
//		ofstream stream;
//		stream.open(file_name);
//		for (size_t i = 0; i < infos.size(); i++) {
//			stream << infos[i] << endl;
//		}
//		stream.close();
//	}
//
//	string ConstructEdgePairFileName(const string output_folder,
//			const string &name, const string &modifier, size_t index) {
//		stringstream ss;
//		ss.clear();
//		ss << output_folder << "/" << name << "_" << index << "_" << modifier
//				<< ".inf";
//		return ss.str();
//	}

//	void WriteWorstEdgesStat(const string &output_folder, double bound) {
//		size_t count = 0;
//		WriteFalseNegativeGaps(output_folder + "/gaps.inf");
//		for (auto iterator = false_positive_infos_.begin();
//				iterator != false_positive_infos_.end(); ++iterator) {
//			if (iterator->weight > bound) {
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp",
//								"histogram", count),
//						pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp",
//								"estimated", count),
//						estimated_pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp", "etalon",
//								count),
//						etalon_pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				count++;
//			}
//		}
//		for (auto iterator = false_negative_infos_.begin();
//				iterator != false_negative_infos_.end(); ++iterator) {
//			if (iterator->weight > bound) {
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp",
//								"histogram", count),
//						pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp",
//								"estimated", count),
//						estimated_pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				WriteEdgePairInfo(
//						ConstructEdgePairFileName(output_folder, "fp", "etalon",
//								count),
//						etalon_pair_info_.GetEdgePairInfo(iterator->first,
//								iterator->second));
//				count++;
//			}
//		}
//	}

};

template<class Graph>
class ClusterStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef PairInfo<EdgeId> Info;
	typedef typename PairedInfoIndex<Graph>::PairInfos Infos;

	const PairedInfoIndex<Graph>& estimated_pair_info_;
	vector<pair<double, double>> weight_variance_stat_;DECL_LOGGER("EstimatedClusterStat")
	;
public:
	ClusterStat(const PairedInfoIndex<Graph>& estimated_pair_info) :
			estimated_pair_info_(estimated_pair_info) {
	}

	virtual ~ClusterStat() {
	}

	virtual void Count() {
		for (auto it = estimated_pair_info_.begin();
				it != estimated_pair_info_.end(); ++it) {
			Infos infos = *it;
			for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
				Info info = *it2;
				//				cerr << "HERE0" <<endl;
				if (gr(info.variance, 0.)) {
					//					cerr << "HERE1" <<endl;
					weight_variance_stat_.push_back(
							make_pair(info.weight, info.variance));
				}
			}
			//todo talk with Anton!!!
			//			for (auto it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
			//				Info info = *it2;
			////				if (gr(info.variance, 0)) {
			//					weight_variance_stat_.push_back(make_pair(info.weight, info.variance));
			////				}
			//			}
		}
		stringstream ss;
		/*	for (auto it = weight_variance_stat_.begin(); it != weight_variance_stat_.end(); ++it) {
		 ss <<*it;//<< "(" << (*it).first << ", " << (*it).second << ")" << " ; ";
		 }
		 */
		copy(weight_variance_stat_.begin(), weight_variance_stat_.end(),
				ostream_iterator<pair<double, double>>(ss, ", "));
		INFO("Estimated cluster stat: " << ss.str());
	}

	vector<pair<double, double>> weight_variance_stat() {
		sort(weight_variance_stat_.begin(), weight_variance_stat_.end());
		return weight_variance_stat_;
	}

};

template<class Graph>
class CoverageStatistics{
    
private:
    Graph& graph_;
    EdgeQuality<Graph> & edge_qual_;
    
    
    bool DoesSuit(VertexId vertex){
        bool ans = true;
        for (size_t i = 0; ans && i<graph_.OutgoingEdgeCount(vertex); i++) 
            ans = ans & math::gr(edge_qual_.quality(graph_.OutgoingEdges(vertex)[i]), 0.);
        for (size_t i = 0; ans && i<graph_.IncomingEdgeCount(vertex); i++) 
            ans = ans & math::gr(edge_qual_.quality(graph_.IncomingEdges(vertex)[i]), 0.);
        return ans; 
    }
    
public:
    CoverageStatistics(Graph& graph, EdgeQuality<Graph>& edge_qual):
    graph_(graph), edge_qual_(edge_qual){
    }

    virtual ~CoverageStatistics(){}

    virtual void Count(){

        map<double, size_t> cov_map;
        map<double, size_t> ratio_map;
        map<double, size_t> len_map;
        size_t area = 0;
        size_t area15 = 0;
        size_t area10 = 0;
        size_t area5 = 0;
        size_t area2 = 0;
        for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter){
            len_map[graph_.length(*iter)]++;
        }
        for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) 
            if (true || DoesSuit(*iter) ){

                double plus_cov = 0.;
                double min_cov = 0.;
                double plus_all_cov = 0.;
                double min_all_cov = 0.;
                bool suit_us = true;

                if (graph_.IncomingEdgeCount(*iter)*graph_.OutgoingEdgeCount(*iter) == 0) continue;
                
                for (size_t i = 0; suit_us && i<graph_.IncomingEdgeCount(*iter); i++) 
                    if (graph_.length(graph_.IncomingEdges(*iter)[i]) < 80){
                        if (math::ge(edge_qual_.quality(graph_.IncomingEdges(*iter)[i]), 1.))
                            plus_cov += graph_.coverage(graph_.IncomingEdges(*iter)[i]);
                        plus_all_cov += graph_.coverage(graph_.IncomingEdges(*iter)[i]);
                    }else suit_us = false;
                for (size_t i = 0; suit_us && i<graph_.OutgoingEdgeCount(*iter); i++) 
                    if (graph_.length(graph_.OutgoingEdges(*iter)[i]) < 80){
                        if (math::ge(edge_qual_.quality(graph_.OutgoingEdges(*iter)[i]), 1.))
                            min_cov += graph_.coverage(graph_.OutgoingEdges(*iter)[i]);
                        min_all_cov += graph_.coverage(graph_.OutgoingEdges(*iter)[i]);
                    }else suit_us = false;

                if (!suit_us) continue;

                if (math::eq(min_cov, 0.) || math::eq(plus_cov, 0.)) continue;

                double delta_cov = math::round(1000.*(plus_cov - min_cov)/(plus_cov + min_cov));
                
                double ratio_cov = math::round(1000.*(plus_cov + min_cov)/(plus_all_cov + min_all_cov));

                if (math::ls(abs(delta_cov), 150.)) area15++;
                if (math::ls(abs(delta_cov), 100.)) area10++;
                if (math::ls(abs(delta_cov), 50.)) area5++;
                if (math::ls(abs(delta_cov), 20.)) area2++;
                area++;
                
                cov_map[delta_cov/10.]++;
                ratio_map[ratio_cov/10.]++;
            
        }

        for (auto iter = ratio_map.begin(); iter != ratio_map.end(); ++iter){
            INFO("Ratio " << (*iter).first << " " << (*iter).second);
        }

        for (auto iter = cov_map.begin(); iter != cov_map.end(); ++iter){
            INFO("Cov " << (*iter).first << " " << (*iter).second);
        }

        INFO("stats_cov "  << area << " " << area2 << " " << area5 << " " << area10 << " " << area15);
        
        for (auto iter = len_map.begin(); iter != len_map.end(); ++iter){
            INFO("Len " << (*iter).first << " " << (*iter).second);
        }
           
    }
    
};

}

#endif /* STATISTICS_HPP_ */

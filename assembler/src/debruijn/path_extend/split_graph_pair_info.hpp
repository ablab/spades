/*
 * split_graph_pair_info.hpp
 *
 *  Created on: May 14, 2013
 *      Author: ira
 */

#ifndef SPLIT_GRAPH_PAIR_INFO_HPP_
#define SPLIT_GRAPH_PAIR_INFO_HPP_

#include "graphio.hpp"
#include "single_threshold_finder.hpp"
using namespace debruijn_graph;

namespace path_extend {
class Basket{

private:
	EdgeId edgeId_;
	size_t index_;

public:
	Basket(EdgeId edgeId, size_t index): edgeId_(edgeId), index_(index){

	}

	Basket(const Basket& b): edgeId_(b.edgeId_), index_(b.index_){

	}

	const EdgeId edgeId() const{
		return edgeId_;
	}

	size_t index() const{
		return index_;
	}

	bool operator<(const Basket& rhs) const {
		if (edgeId() != rhs.edgeId()) {
			return edgeId() < rhs.edgeId();
		}
		return index() < rhs.index();
	}

	bool operator==(const Basket& rhs) const {
		return edgeId() == rhs.edgeId() && index() == rhs.index();
	}
};

struct PairInfo{
	double weight_;
	double distance_;
	size_t count_;

	PairInfo(){
		weight_ = 0.;
		distance_ = 0.;
		count_ = 0;
	}

	PairInfo(double weight, double distance, size_t count = 0): weight_(weight), distance_(distance), count_(count){

	}

};

class EdgePairInfo{

private:
	EdgeId edgeId_;
	size_t basket_size_;
	vector<map<Basket, PairInfo> > pair_info_;

public:
	EdgePairInfo(){
		basket_size_ = 0;
	}

	EdgePairInfo(conj_graph_pack& gp, EdgeId edgeId, size_t basket_size):
		edgeId_(edgeId), basket_size_(basket_size){
		size_t count_baskets = gp.g.length(edgeId) / basket_size_ + 1;
		for (size_t index = 0; index < count_baskets; ++index){
			pair_info_.push_back(map<Basket, PairInfo>());
		}
	}

	EdgePairInfo(const EdgePairInfo& pairInfo):
		edgeId_(pairInfo.edgeId_), basket_size_(pairInfo.basket_size_){
			for (size_t index = 0; index < pairInfo.pair_info_.size(); ++index){
				pair_info_.push_back(pairInfo.pair_info_[index]);
			}
	}

	void AddPairInfo(size_t pos_begin1, size_t pos_end1,
			EdgeId edgeId2, size_t pos_begin2, size_t pos_end2,
			double weight, double edge_distance) {
		size_t begin_basket_index1 = GetBasketIndex(pos_begin1);
		size_t end_basket_index1 = GetBasketIndex(pos_end1);
		size_t begin_basket_index2 = GetBasketIndex(pos_begin2);
		size_t end_basket_index2 = GetBasketIndex(pos_end2);
		for (size_t index1 = begin_basket_index1; index1 <= end_basket_index1; ++index1){
			for (size_t index2 = begin_basket_index2; index2 <= end_basket_index2; ++index2){
				AddPairInfoToBasket(index1, edgeId2, index2, weight, edge_distance);
			}
		}
	}

	map<Basket, PairInfo>& GetBasketInfo(size_t index){
		return pair_info_.at(index);
	}

	size_t size(){
		return pair_info_.size();
	}

private:
	size_t GetBasketIndex(size_t pos) const{
		return pos / basket_size_;
	}

	void AddPairInfoToBasket(size_t index1,
			EdgeId edgeId2, size_t index2, double weight, double edge_distance) {
		Basket basket2(edgeId2, index2);
		if (pair_info_[index1].find(basket2) == pair_info_[index1].end()) {
			pair_info_[index1][basket2] = PairInfo(0.0, 0);
		}
		PairInfo oldPairInfo = pair_info_[index1][basket2];
		double basket_distance = GetBasketDistance(edge_distance, index1, index2);
		pair_info_[index1][basket2] = PairInfo(oldPairInfo.weight_ + weight,
				CountNewDistance(oldPairInfo, basket_distance),
				oldPairInfo.count_ + 1);
	}

	double CountNewDistance(PairInfo& oldPairInfo, double distance){
		return (oldPairInfo.distance_ * oldPairInfo.count_ + distance) / (oldPairInfo.count_ + 1);
	}

	double GetBasketDistance(double edge_distance, double index1, double index2){
		return edge_distance - index1 * basket_size_ + index2 * basket_size_;
	}
};

class BasketsPairInfoIndex{

private:
	conj_graph_pack& gp_;
	size_t basket_size_;
	map<EdgeId, EdgePairInfo> pair_info_;

public:
	BasketsPairInfoIndex(conj_graph_pack& gp, size_t basket_size): gp_(gp), basket_size_(basket_size){

	}

	void AddPairInfo(EdgeId edgeId1, size_t pos_begin1, size_t pos_end1,
			EdgeId edgeId2, size_t pos_begin2, size_t pos_end2,
			double weight, double edge_distance)  {
		if (pair_info_.find(edgeId1) == pair_info_.end()){
			EdgePairInfo edgePairInfo2(gp_, edgeId1, basket_size_);
			pair_info_.insert(make_pair(edgeId1, edgePairInfo2));
		}
		pair_info_[edgeId1].AddPairInfo(pos_begin1, pos_end1, edgeId2, pos_begin2, pos_end2, weight, edge_distance);
	}

	EdgePairInfo& GetEdgePairInfo(EdgeId edgeId){
		return pair_info_[edgeId];
	}

};

class SplitGraphPairInfo{

private:
	conj_graph_pack& gp_;
	PairedInfoLibrary& lib_;
	size_t lib_index_;
	size_t basket_size_;
	BasketsPairInfoIndex basketIndex_;
public:
	SplitGraphPairInfo(conj_graph_pack& gp, PairedInfoLibrary& lib, size_t lib_index, size_t basket_size):
		gp_(gp), lib_(lib), lib_index_(lib_index), basket_size_(basket_size), basketIndex_(gp, basket_size){

	}

	void ProcessReadPairs() {
		auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(cfg::get().ds.reads[lib_index_], true, cfg::get().ds.reads[lib_index_].data().mean_insert_size);
		SingleStreamType paired_streams(paired_stream.get());
		io::IReader<io::PairedRead>& stream = paired_streams.back();
		stream.reset();

		while (!stream.eof()) {
		    io::IReader<io::PairedRead>::read_type p_r;
		    stream >> p_r;
		    Sequence read1 = p_r.first().sequence();
		    Sequence read2 = p_r.second().sequence();
			size_t read_distance = p_r.distance();
			NewExtendedSequenceMapper<Graph> mapper(gp_.g, gp_.index,
					gp_.kmer_mapper, gp_.g.k() + 1);
			MappingPath<EdgeId> path1 = mapper.MapSequence(read1);
			MappingPath<EdgeId> path2 = mapper.MapSequence(read2);

			for (size_t i = 0; i < path1.size(); ++i) {
				pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];

				for (size_t j = 0; j < path2.size(); ++j) {
					pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
					double weight = PairedReadCountWeight(mapping_edge_1.second,
							mapping_edge_2.second);

					size_t kmer_distance = read_distance
							+ mapping_edge_2.second.initial_range.end_pos
							- mapping_edge_1.second.initial_range.start_pos;
					int edge_distance = kmer_distance
							+ mapping_edge_1.second.mapped_range.start_pos
							- mapping_edge_2.second.mapped_range.end_pos;

					basketIndex_.AddPairInfo(mapping_edge_1.first,
							mapping_edge_1.second.mapped_range.start_pos,
							mapping_edge_1.second.mapped_range.end_pos,
							mapping_edge_2.first,
							mapping_edge_2.second.mapped_range.start_pos,
							mapping_edge_2.second.mapped_range.end_pos,
							weight, (double)edge_distance
							);
				}
			}
			//INFO(++n);
		}
	}

	double FindThreshold(double min_length_long_edge, int insert_size_min, int insert_size_max) {
		Graph& graph = gp_.g;
		vector<double> good_pi;
		vector<double> bad_pi;

		for (auto edge = graph.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
			if (graph.length(*edge) > min_length_long_edge) {
				if (graph.int_id(*edge) <= 0) {
					continue;
				}
				EdgePairInfo& edgePairInfo = basketIndex_.GetEdgePairInfo(*edge);
				if (edgePairInfo.size() == 0){
					continue;
				}

				size_t countBackets = LastTestBasketIndex(*edge, insert_size_max, edgePairInfo);
				for (size_t index = 0; index <= countBackets; ++index){
					map<Basket, PairInfo>& basketInfo = edgePairInfo.GetBasketInfo(index);
					set<size_t> pairBaskets = GetPairBaskets(index, insert_size_min, insert_size_max, edgePairInfo);

					for (auto pairBasketIter = basketInfo.begin(); pairBasketIter != basketInfo.end(); ++pairBasketIter){
						PairInfo& pairInfo = pairBasketIter->second;
						if (pairBasketIter->first.edgeId() == *edge && pairBaskets.find(pairBasketIter->first.index()) != pairBaskets.end()){
							good_pi.push_back(GetNormalizedWeight(pairInfo));
						} else {
							bad_pi.push_back(GetNormalizedWeight(pairInfo));
						}
					}

				}
			}
		}
		double threshold = find_intersection(good_pi, bad_pi);
		INFO("WE FOUND THRESHOLD " << threshold <<" good_pi_size " << good_pi.size() << " bad_pi_size " << bad_pi.size());
		return threshold;
	}

private:

	size_t LastTestBasketIndex(EdgeId edgeId, int insert_size_max, EdgePairInfo& edgePairInfo){
		return min((gp_.g.length(edgeId) - insert_size_max) / basket_size_, edgePairInfo.size() - 1);
	}

	size_t FindBeginPairBasket(size_t index, int insert_size_min, EdgePairInfo& edgePairInfo){
		return min(index + insert_size_min / basket_size_, edgePairInfo.size() - 1);
	}

	size_t FindEndPairBasket(size_t index, int insert_size_max, EdgePairInfo& edgePairInfo){
		return min(index + insert_size_max/ basket_size_, edgePairInfo.size() - 1);
	}

	set<size_t> GetPairBaskets(size_t index, int insert_size_min, int insert_size_max, EdgePairInfo& edgePairInfo){
		set<size_t> result;
		size_t begin = FindBeginPairBasket(index, insert_size_min, edgePairInfo);
		size_t end = FindEndPairBasket(index, insert_size_max, edgePairInfo);
		for (size_t pair_index = begin; pair_index <= end; ++pair_index){
			result.insert(pair_index);
		}
		return result;
	}

	double GetNormalizedWeight(PairInfo& pi){
		//INFO("weight "<< pi.weight_ << " distance " << pi.distance_ << " " << pi.weight_ / lib_.IdealPairedInfo(basket_size_, basket_size_, pi.distance_));
		return pi.weight_ / lib_.IdealPairedInfo(basket_size_, basket_size_, pi.distance_);
	}
};

} /* path_extend */

#endif /* SPLIT_GRAPH_PAIR_INFO_HPP_ */

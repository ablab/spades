#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "sequence.hpp"
#include <cmath>
#include <map>

#define MERGE_DATA_ABSOLUTE_DIFFERENCE 0
//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace de_bruijn {

template<class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const int max_difference_;

public:

	class PairInfo {
		friend class PairedInfoIndex;
		friend class PairedInfoIndexData;
		friend class EdgePairIterator;
	private:
		EdgeId first_;
		EdgeId second_;
		double d_;//distance between starts. Can be negative
		double weight_;
	public:
		EdgeId first() const {
			return first_;
		}
		EdgeId second() const {
			return second_;
		}

		int d() const {
			int res = (int) (std::abs(d_) + 0.5 + 1e-9);
			if (d_ < 0)
				res = -res;
			return res;
		}

		int exact_d() const {
			return d_;
		}

		double weight() const {
			return weight_;
		}

		PairInfo(EdgeId first, EdgeId second, double d, double weight) :
			first_(first), second_(second), d_(d), weight_(weight) {
		}

		bool operator==(const PairInfo& rhs) const {
			return first_ == rhs.first_ && second_ == rhs.second_ && d_
					== rhs.d_/* && weight_ == rhs.weight_*/;
		}
	};

	typedef vector<PairInfo> PairInfos;
	typedef typename PairInfos::const_iterator infos_iterator;

private:
	//todo try storing set<PairInfo>
	class PairInfoIndexData {
	public:
		typedef multimap<pair<EdgeId, EdgeId> , pair<double, double>> Data;
		typedef typename Data::iterator data_iterator;
		typedef typename Data::const_iterator const_data_iterator;
	private:
		Data data_;

		PairInfo BackwardInfo(const PairInfo& pair_info) {
			return PairInfo(pair_info.second_, pair_info.first_, -pair_info.d_,
					pair_info.weight_);
		}

		const PairInfo AsPairInfo(
				const pair<pair<EdgeId, EdgeId> , pair<double, double>>& pair) {
			return PairInfo(pair.first.first, pair.first.second,
					pair.second.first, pair.second.second);
		}

		const pair<EdgeId, EdgeId> EdgePair(const PairInfo& pair_info) {
			return make_pair(pair_info.first_, pair_info.second_);
		}

		const pair<pair<EdgeId, EdgeId> , pair<double, double>> AsPairOfPairs(
				const PairInfo& pair_info) {
			return make_pair(EdgePair(pair_info),
					make_pair(pair_info.d_, pair_info.weight_));
		}

		void UpdateSingleInfo(const PairInfo& info, const double d,
				const double weight) {
			bool updated = false;
			for (data_iterator lower = data_.lower_bound(EdgePair(info)),
					upper = data_.upper_bound(EdgePair(info)); lower != upper; ++lower) {
				const PairInfo& existing_info = AsPairInfo(*lower);
				if (existing_info == info) {
					lower->second.first = d;
					lower->second.second = weight;
					updated = true;
					break;
				}
			}
			assert(updated);
		}

	public:
		data_iterator begin() {
			return data_.begin();
		}

		data_iterator end() {
			return data_.end();
		}

		void AddPairInfo(const PairInfo& pair_info) {
			if (pair_info.first_ == pair_info.second_ && pair_info.d_ == 0) {
				data_.insert(AsPairOfPairs(pair_info));
			} else {
				data_.insert(AsPairOfPairs(pair_info));
				data_.insert(AsPairOfPairs(BackwardInfo(pair_info)));
			}
		}

		void DeleteEdgeInfo(EdgeId e) {
			set<EdgeId> paired_edges;
			for (const_data_iterator lower = LowerBound(e), upper = UpperBound(
					e); lower != upper; ++lower) {
				paired_edges.insert((*lower).first.second);
			}
			for (typename set<EdgeId>::const_iterator it = paired_edges.begin(); it
					!= paired_edges.end(); ++it) {
				data_.erase(make_pair(*it, e));
			}
			data_.erase(LowerBound(e), UpperBound(e));
		}

		PairInfos GetEdgeInfos(EdgeId e) {
			vector<PairInfo> answer;
			for (const_data_iterator lower = LowerBound(e), upper = UpperBound(
					e); lower != upper; ++lower) {
				answer.push_back(AsPairInfo(*lower));
			}
			return answer;
		}

		PairInfos GetEdgePairInfos(EdgeId e1, EdgeId e2) {
			vector<PairInfo> answer;
			for (const_data_iterator lower = data_.lower_bound(
					make_pair(e1, e2)), upper = data_.upper_bound(
					make_pair(e1, e2)); lower != upper; ++lower) {
				answer.push_back(AsPairInfo(*lower));
			}
			return answer;
		}

		void UpdateInfo(const PairInfo& info, const int d, const double weight) {
			if (info.first() == info.second() && d == 0) {
				UpdateSingleInfo(info, d, weight);
			} else {
				UpdateSingleInfo(info, d, weight);
				UpdateSingleInfo(BackwardInfo(info), -d, weight);
			}
		}

		void clear() {
			data_.clear();
		}

		data_iterator LowerBound(EdgeId e) {
			return data_.lower_bound(make_pair(e, (EdgeId) 0));
		}

		data_iterator UpperBound(EdgeId e) {
			return data_.upper_bound(make_pair(e, (EdgeId) ((size_t) -1)));
		}

		data_iterator LowerBound(EdgeId e1, EdgeId e2) {
			return data_.lower_bound(make_pair(e1, e2));
		}

		data_iterator UpperBound(EdgeId e1, EdgeId e2) {
			return data_.upper_bound(make_pair(e1, e2));
		}
	};

public:
	class EdgePairIterator {
		typename PairInfoIndexData::data_iterator position_;
		PairedInfoIndex<Graph> &index_;
	public:
		EdgePairIterator(typename PairInfoIndexData::data_iterator position,
				PairedInfoIndex<Graph> &index) :
			position_(position), index_(index) {
		}

		bool operator==(const EdgePairIterator &other) {
			return this->position_ == other.position_;
		}

		bool operator!=(const EdgePairIterator &other) {
			return this->position_ != other.position_;
		}

		PairInfos operator*() const {
			pair<EdgeId, EdgeId> currentPair = position_->first;
			return index_.GetEdgePairInfo(currentPair.first, currentPair.second);
		}

		void operator++() {
			pair<EdgeId, EdgeId> currentPair = position_->first;
			position_ = index_.data_.UpperBound(currentPair.first,
					currentPair.second);
		}
	};

	EdgePairIterator begin() {
		return EdgePairIterator(data_.begin(), *this);
	}

	EdgePairIterator end() {
		return EdgePairIterator(data_.end(), *this);
	}

	//begin-end insert size supposed
	PairedInfoIndex(Graph &g,
			int max_difference = MERGE_DATA_ABSOLUTE_DIFFERENCE) :
		GraphActionHandler<Graph> ("PairedInfoIndex"),
				max_difference_(max_difference), graph_(g) {
		g.AddActionHandler(this);
	}

	virtual ~PairedInfoIndex() {
		graph_.RemoveActionHandler(this);
	}

	template<size_t kmer_size, class Stream>
	void FillIndex(const EdgeIndex<kmer_size + 1, Graph>& index, Stream& stream) {
		data_.clear();
		auto it = graph_.SmartEdgeBegin();
		for (auto it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it; ++it) {
			AddPairInfo(PairInfo(*it, *it, 0, 1));
		}
		typedef Seq<kmer_size + 1> KPOMer;
		de_bruijn::SimpleSequenceMapper<kmer_size, Graph> read_threader(graph_,
				index);
		while (!stream.eof()) {
			PairedRead p_r;
			stream >> p_r;
			ProcessPairedRead(p_r, read_threader);
		}
	}

	double sum() {
		double res = 0;
		for (auto it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it; ++it)
			for (auto it1 = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd()
					!= it1; ++it1) {
				PairInfos vec = GetEdgePairInfo(*it1, *it);
				if (vec.size() != 0) {
					for (size_t i = 0; i < vec.size(); i++)
						res += vec[i].weight_;
				}
			}
		return res;
	}

private:

	template<size_t kmer_size>
	void ProcessPairedRead(const PairedRead& p_r,
			de_bruijn::SimpleSequenceMapper<kmer_size, Graph> &read_threader) {
		Sequence read1 = p_r.first().getSequence();
		Sequence read2 = p_r.second().getSequence();
		de_bruijn::Path<EdgeId> path1 = read_threader.MapSequence(read1);
		de_bruijn::Path<EdgeId> path2 = read_threader.MapSequence(read2);
		size_t distance = CountDistance(p_r);
		int current_distance1 = distance + path1.start_pos()
				- path2.start_pos();
		for (size_t i = 0; i < path1.size(); ++i) {
			int current_distance2 = current_distance1;
			for (size_t j = 0; j < path2.size(); ++j) {
				//				double weight = CorrectLength(path1, i) * CorrectLength(path2,
				//						j);
				double weight = 1;
				PairInfo
						new_info(path1[i], path2[j], current_distance2, weight);
				AddPairInfo(new_info);
				current_distance2 += graph_.length(path2[j]);
			}
			current_distance1 -= graph_.length(path1[i]);
		}
	}

	Graph& graph_;
	PairInfoIndexData data_;

	size_t CorrectLength(const de_bruijn::Path<EdgeId>& path, size_t index) {
		size_t result = graph_.length(path[index]);
		if (index == 0) {
			result -= path.start_pos();
		}
		if (index == path.size() - 1) {
			result -= graph_.length(path[index]) - path.end_pos();
		}
		return result;
	}

	void PassEdge(size_t edge_length, size_t &path_nucls_passed) {
		if (path_nucls_passed == 0) {
			path_nucls_passed += graph_.k();
		}
		path_nucls_passed += edge_length;
	}

	size_t CountDistance(const PairedRead& paired_read) {
		return paired_read.distance() - paired_read.second().size();
	}

	bool CanMergeData(const PairInfo& info1, const PairInfo& info2) {
		if (info1.first_ != info2.first_ || info1.second_ != info2.second_)
			return false;
		if (std::abs(info2.d_ - info1.d_) <= max_difference_) {
			return true;
		}
		return false;
	}

	void MergeData(const PairInfo& info1, const PairInfo& info2) {
		assert(info1.first_ == info2.first_ && info1.second_ == info2.second_);
		double newWeight = info1.weight_ + info2.weight_;
		double newD = (info1.d_ * info1.weight_ + info2.d_ * info2.weight_)
				/ newWeight;
		if (info1.first_ == info1.second_ && (info1.d_ == 0 || info2.d_ == 0)) {
			newD = 0;
			newWeight += info2.weight_;
		}
		data_.UpdateInfo(info1, newD, newWeight);
	}

public:
	void AddPairInfo(const PairInfo& pair_info) {
		PairInfos pair_infos = data_.GetEdgePairInfos(pair_info.first_,
				pair_info.second_);
		for (infos_iterator it = pair_infos.begin(); it != pair_infos.end(); ++it) {
			if (CanMergeData(*it, pair_info)) {
				MergeData(*it, pair_info);
				return;
			}
		}
		data_.AddPairInfo(pair_info);
	}

private:
	void RemoveEdgeInfo(EdgeId edge) {
		data_.DeleteEdgeInfo(edge);
	}

	void OutputEdgeData(EdgeId edge1, EdgeId edge2, ostream &os = cout) {
		PairInfos vec = GetEdgePairInfo(edge1, edge2);
		if (vec.size() != 0) {
			os << edge1 << " " << graph_.length(edge1) << " " << edge2 << " "
					<< graph_.length(edge2) << endl;
			if (graph_.EdgeEnd(edge1) == graph_.EdgeStart(edge2))
				os << "+" << endl;
			if (graph_.EdgeEnd(edge2) == graph_.EdgeStart(edge1))
				os << "-" << endl;
			int min = INT_MIN;
			for (size_t i = 0; i < vec.size(); i++) {
				int next = -1;
				for (size_t j = 0; j < vec.size(); j++) {
					if (vec[j].d() > min && (next == -1 || vec[next].d()
							> vec[j].d())) {
						next = j;
					}
				}
				os << vec[next].d() << " " << vec[next].weight() << endl;
				if (next == -1) {
					assert(false);
				}
				if (vec[next].d() > 100000) {
					assert(false);
				}
				min = vec[next].d();
			}
		}
	}

	void TransferInfo(EdgeId old_edge, EdgeId new_edge, int shift = 0,
			double weight_scale = 1.0) {
		PairInfos pair_infos = GetEdgeInfo(old_edge);
		for (size_t j = 0; j < pair_infos.size(); ++j) {
			PairInfo old_pair_info = pair_infos[j];
			if (old_edge != old_pair_info.second()) {
				AddPairInfo(
						PairInfo(new_edge, old_pair_info.second_,
								old_pair_info.d_ - shift,
								weight_scale * old_pair_info.weight_));
			} else {
				AddPairInfo(
						PairInfo(new_edge, new_edge, old_pair_info.d_,
								weight_scale * 0.5 * old_pair_info.weight_));
			}
		}
		RemoveEdgeInfo(old_edge);
	}

public:

	void OutputData(ostream &os = cout) {
		for (auto it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it; ++it)
			for (auto it1 = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd()
					!= it1; ++it1) {
				OutputEdgeData(*it, *it1, os);
			}
	}
	void OutputData(string fileName) {
		ofstream s;
		s.open(fileName.c_str());
		OutputData(s);
		s.close();
	}

	PairInfos GetEdgeInfo(EdgeId edge) {
		return data_.GetEdgeInfos(edge);
	}

	PairInfos GetEdgePairInfo(EdgeId first, EdgeId second) {
		return data_.GetEdgePairInfos(first, second);
	}

	virtual void HandleAdd(EdgeId e) {
		this->AddPairInfo(PairInfo(e, e, 0, 1));
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}

	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
		this->AddPairInfo(PairInfo(new_edge, new_edge, 0, 1));
		int shift = 0;
		for (size_t i = 0; i < old_edges.size(); ++i) {
			EdgeId old_edge = old_edges[i];
			TransferInfo(old_edge, new_edge, shift);
			shift -= graph_.length(old_edge);
		}
	}

	virtual void HandleGlue(EdgeId old_edge, EdgeId new_edge) {
		TransferInfo(old_edge, new_edge);
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		double prop = (double) graph_.length(new_edge1) / graph_.length(
				old_edge);
		size_t shift = 0;
		TransferInfo(old_edge, new_edge1, shift, prop);
		PassEdge(graph_.length(new_edge1), shift);
		TransferInfo(old_edge, new_edge1, shift, 1 - prop);
	}

};

template<class Graph>
class SimpleOfflineClusterer {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename PairedInfoIndex<Graph>::PairInfo PairInfo;
	typedef vector<PairInfo> PairInfos;

public:
	PairedInfoIndex<Graph> &not_clustered_;

	SimpleOfflineClusterer(PairedInfoIndex<Graph> &not_clustered) :
		not_clustered_(not_clustered) {
	}

	typename PairedInfoIndex<Graph>::PairInfos ProcessEdgePair(
			const typename PairedInfoIndex<Graph>::PairInfos &infos) {
		EdgeId edge1 = infos[0].first();
		EdgeId edge2 = infos[0].second();
		double d_sum = 0;
		double weight_sum = 0;
		for (size_t i = 0; i < infos.size(); i++) {
			d_sum += infos[i].exact_d();
			weight_sum += infos[i].weight();
		}
		PairInfo sum_info(edge1, edge2, d_sum / infos.size(), weight_sum);
		PairInfos result;
		result.push_back(sum_info);
		return result;
	}

	void cluster(PairedInfoIndex<Graph> &clustered) {
		assert(&not_clustered_ != &clustered);
		for (typename PairedInfoIndex<Graph>::EdgePairIterator it =
				not_clustered_.begin(); it != not_clustered_.end(); ++it) {
			typename PairedInfoIndex<Graph>::PairInfos newInfos =
					ProcessEdgePair(*it);
			for (size_t i = 0; i < newInfos.size(); i++) {
				clustered.AddPairInfo(newInfos[i]);
			}
		}
	}
};

}

#endif /* PAIRED_INFO_HPP_ */

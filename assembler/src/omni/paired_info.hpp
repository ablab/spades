#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
//#include "utils.hpp"
//#include "sequence.hpp"
//#include "paired_read.hpp"
#include "common/io/paired_read.hpp"
#include <cmath>
#include <map>
#include <limits>

//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3
#define E 1e-6
namespace omnigraph {

/**
 * PairInfo class represents basic data unit for paired information: edges first and second appear
 * in genome at distance d_ and this information has weight weight_.
 */
template<typename EdgeId>
struct PairInfo {
	const EdgeId first;
	const EdgeId second;
	const double d; //distance between starts. Can be negative
	const double weight;
	const double variance;

	PairInfo(const PairInfo& pair_info) :
		first(pair_info.first), second(pair_info.second), d(pair_info.d),
				weight(pair_info.weight), variance(pair_info.variance) {

	}

	PairInfo(EdgeId first, EdgeId second, double d, double weight,
			double variance = 0.) :
		first(first), second(second), d(d), weight(weight), variance(variance) {
	}

	const PairInfo set_first(EdgeId first) const {
		return PairInfo(first, this->second, this->d, this->weight,
				this->variance);
	}

	const PairInfo set_second(EdgeId second) const {
		return PairInfo(this->first, second, this->d, this->weight,
				this->variance);
	}

	const PairInfo set_distance(double d) const {
		return PairInfo(this->first, this->second, d, this->weight,
				this->variance);
	}

	const PairInfo set_weight(EdgeId first) const {
		return PairInfo(this->first, this->second, this->d, weight,
				this->variance);
	}

	const PairInfo set_variance(EdgeId first) const {
		return PairInfo(this->first, this->second, this->d, this->weight,
				variance);
	}

	bool operator<(const PairInfo& rhs) const {
		const PairInfo &lhs = *this;
		return lhs.first == rhs.first ? lhs.second == rhs.second ? lhs.d + E
				< rhs.d : lhs.second < rhs.second : lhs.first < rhs.first;
	}

	const PairInfo& operator=(const PairInfo& pair_info) {
		if (this != &pair_info) {
			assert(false);
		}
		return *this;
	}

	/**
	 * Two paired infos are considered equal if they coinside in all parameters except for the weight of
	 * info.
	 */
	bool operator==(const PairInfo& rhs) const {
		const PairInfo &lhs = *this;
		return !(lhs < rhs || rhs < lhs);

		// uncomment this for speed-up:

		//    return lhs.first  == rhs.first      &&
		//           lhs.second == rhs.second     &&
		//           lhs.d      == rhs.d     /*   &&
		//           lhs.weight == rhs.weight*/;
	}
};

//typedef vector<PairInfo<> > PairInfos;

template<typename EdgeId>
PairInfo<EdgeId> MinPairInfo(EdgeId id) {
	return PairInfo<EdgeId> (id, (EdgeId) 0/*numeric_limits<EdgeId>::min()*/,
			-100000000/*numeric_limits<double>::min()*/, 0.);
}

template<typename EdgeId>
PairInfo<EdgeId> MaxPairInfo(EdgeId id) {
	return PairInfo<EdgeId> (id,
			(EdgeId) (-1)/*numeric_limits<EdgeId>::max()*/,
			1000000000/*numeric_limits<double>::max()*/, 0.);
}

template<typename EdgeId>
PairInfo<EdgeId> MinPairInfo(EdgeId e1, EdgeId e2) {
	return MinPairInfo(e1).set_second(e2);
}

template<typename EdgeId>
PairInfo<EdgeId> MaxPairInfo(EdgeId e1, EdgeId e2) {
	return MaxPairInfo(e1).set_second(e2);
}

/**
 * Method returns approximate distance between occurrences of edges in genome rounded to the nearest
 * integer. In case of a tie closest to 0 value is chosen thus one can assume that distance
 * is rounded the same way as opposite one.
 */
template<typename EdgeId>
int rounded_d(PairInfo<EdgeId> const& pi) {
	int res = (int) (std::abs(pi.d) + 0.5 + 1e-9);
	if (pi.d < 0)
		res = -res;
	return res;
}

template<typename EdgeId>
PairInfo<EdgeId> BackwardInfo(const PairInfo<EdgeId>& pi) {
	return PairInfo<EdgeId> (pi.second, pi.first, -pi.d, pi.weight);
}

template<typename EdgeId>
bool IsSymmetric(PairInfo<EdgeId> const& pi) {
	return pi.first == pi.second && pi.d == 0;
}

template<typename EdgeId>
void PrintPairInfo(const PairInfo<EdgeId>& pair_info) {
	cerr << pair_info.first << "  " << pair_info.second << "  " << pair_info.d
			<< "  " << pair_info.weight << endl;
}

template<typename EdgeId>
void PrintPairInfos(const vector<PairInfo<EdgeId>>& pair_infos) {
	for (auto it = pair_infos.begin(); it != pair_infos.end(); ++it) {
		PrintPairInfo(*it);
	}
}

//todo try storing set<PairInfo>
template<typename EdgeId>
class PairInfoIndexData {
public:
	typedef set<PairInfo<EdgeId>> Data;
	typedef typename Data::iterator data_iterator;
	typedef typename Data::const_iterator data_const_iterator;
	typedef vector<PairInfo<EdgeId>> PairInfos;

	typedef std::pair<data_const_iterator, data_const_iterator> iterator_range;

public:
	void UpdateSingleInfo(const PairInfo<EdgeId>& info, double d, double weight) {
		size_t count = data_.erase(info);
		assert(count != 0);
		data_.insert(PairInfo<EdgeId> (info.first, info.second, d, weight));
	}

	void ReplaceFirstEdge(const PairInfo<EdgeId>& info, EdgeId newId) {
		//		size_t count = data_.erase(info);
		//	assert(count != 0);
		data_.insert(PairInfo<EdgeId> (newId, info.second, info.d, info.weight));
	}
public:
	data_iterator begin() {
		auto itp = data_.begin();
		int cnt = 1;
		for (auto it = data_.begin(); it != data_.end(); ++it) {
			if (it->first == itp->first && it->second == itp->second && it
					!= data_.begin()) {
				cnt++;
			} else {
				if (it != data_.begin()) {
					cnt = 1;
				}
			}
			itp = it;
		}
		return data_.begin();
	}

	data_iterator end() {
		return data_.end();
	}

	size_t size() {
		return data_.size();
	}

	void AddPairInfo(const PairInfo<EdgeId>& pair_info, bool addSymmetric = 1) {
		TRACE("REALLY ADD:" << pair_info.first << pair_info.second);

		data_.insert(pair_info);

		if (!IsSymmetric(pair_info) && addSymmetric)
			data_.insert(BackwardInfo(pair_info));
	}

	void DeleteEdgeInfo(EdgeId e) {
		set<PairInfo<EdgeId>> paired_edges;

		for (auto lower = LowerBound(e), upper = UpperBound(e); lower != upper; ++lower) {
			paired_edges.insert(BackwardInfo(*lower));
		}

		for (auto it = paired_edges.begin(); it != paired_edges.end(); ++it) {
			data_.erase(*it);
		}

		data_.erase(LowerBound(e), UpperBound(e));
	}

	PairInfos GetEdgeInfos(EdgeId e) {
		return PairInfos(LowerBound(e), UpperBound(e));
	}

	PairInfos GetEdgePairInfos(EdgeId e1, EdgeId e2) {
		return PairInfos(LowerBound(e1, e2), UpperBound(e1, e2));
	}

	void UpdateInfo(const PairInfo<EdgeId>& info, const int d,
			const double weight) {
		UpdateSingleInfo(info, d, weight);

		if (!IsSymmetric(info))
			UpdateSingleInfo(BackwardInfo(info), -d, weight);
	}

	void clear() {
		data_.clear();
	}

	data_iterator LowerBound(EdgeId e) const {
		return data_.lower_bound(MinPairInfo(e));
	}

	data_iterator UpperBound(EdgeId e) {
		return data_.upper_bound(MaxPairInfo(e));
	}

	data_iterator LowerBound(EdgeId e1, EdgeId e2) {
		return data_.lower_bound(MinPairInfo(e1, e2));
	}

	data_iterator UpperBound(EdgeId e1, EdgeId e2) {
		return data_.upper_bound(MaxPairInfo(e1, e2));
	}

private:
	Data data_;
};

/**
 * PairedInfoIndex stores information about edges connected by paired reads and synchronizes this info with
 * graph.
 */
template<class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {

public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<PairInfo<EdgeId>> PairInfos;

private:
	const int max_difference_;

public:
	/**
	 * Class EdgePairIterator is used to iterate through paires of edges which have information about distance
	 * between them stored in PairedInfoIndex.
	 */
	class EdgePairIterator {
		typename PairInfoIndexData<EdgeId>::data_iterator position_;
		PairedInfoIndex<Graph> &index_;
	public:
		EdgePairIterator(
				typename PairInfoIndexData<EdgeId>::data_iterator position,
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
			return index_.GetEdgePairInfo(position_->first, position_->second);
		}

		void operator++() {
			position_ = index_.data_.UpperBound(position_->first,
					position_->second);
		}
	};

	EdgePairIterator begin() {
		return EdgePairIterator(data_.begin(), *this);
	}

	EdgePairIterator end() {
		return EdgePairIterator(data_.end(), *this);
	}

	//begin-end insert size supposed
	PairedInfoIndex(Graph &g, int max_difference = 0) :
		GraphActionHandler<Graph> ("PairedInfoIndex"),
				max_difference_(max_difference), graph_(g) {
		g.AddActionHandler(this);
	}

	virtual ~PairedInfoIndex() {
		TRACE("~PairedInfoIndex");
		graph_.RemoveActionHandler(this);
		TRACE("~PairedInfoIndex ok");
	}

	double sum() {
		double res = 0;
		for (auto it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it; ++it)
			for (auto it1 = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd()
					!= it1; ++it1) {
				PairInfos vec = GetEdgePairInfo(*it1, *it);
				if (vec.size() != 0) {
					for (size_t i = 0; i < vec.size(); i++)
						res += vec[i].weight;
				}
			}
		return res;
	}

private:

	Graph& graph_;
	PairInfoIndexData<EdgeId> data_;

	size_t CorrectLength(const Path<EdgeId>& path, size_t index) {
		size_t result = graph_.length(path[index]);
		if (index == 0) {
			result -= path.start_pos();
		}
		if (index == path.size() - 1) {
			result -= graph_.length(path[index]) - path.end_pos();
		}
		return result;
	}

	//	void PassEdge(size_t edge_length, size_t &path_nucls_passed) {
	//		if (path_nucls_passed == 0) {
	//			path_nucls_passed += graph_.k();
	//		}
	//		path_nucls_passed += edge_length;
	//	}

	bool CanMergeData(const PairInfo<EdgeId>& info1,
			const PairInfo<EdgeId>& info2) {
		if (info1.first != info2.first || info1.second != info2.second)
			return false;
		if (std::abs(info2.d - info1.d) <= max_difference_) {
			return true;
		}
		return false;
	}

	void MergeData(const PairInfo<EdgeId>& info1, const PairInfo<EdgeId>& info2) {
		assert(info1.first == info2.first && info1.second == info2.second);
		double newWeight = info1.weight + info2.weight;
		double newD = (info1.d * info1.weight + info2.d * info2.weight)
				/ newWeight;
		if (info1.first == info1.second && (info1.d == 0 || info2.d == 0)) {
			newD = 0;
			newWeight += info2.weight;
		}
		data_.UpdateInfo(info1, newD, newWeight);
	}

	int NearestClusterIndex(const vector<PairInfo<EdgeId>>& current_pair_infos,
			const PairInfo<EdgeId>& new_info) {
		double min_dist = max_difference_ + 1e-9;
		int answer = -1;
		for (size_t i = 0; i < current_pair_infos.size(); ++i) {
			if (std::abs(new_info.d - current_pair_infos[i].d) < min_dist
					+ 1e-9)
				if (std::abs(new_info.d - current_pair_infos[i].d) < min_dist
						- 1e-9 || answer == -1 || std::abs(
						current_pair_infos[answer].d) > std::abs(
						current_pair_infos[i].d)) {
					min_dist = std::abs(new_info.d - current_pair_infos[i].d);
					answer = i;
				}
		}
		return answer;
	}

public:
	/**
	 * Method allows to add pair info to index directly instead of filling it from stream.
	 */
	void AddPairInfo(const PairInfo<EdgeId>& pair_info, bool add_reversed = 1) {
		TRACE(
				"IN ADD:" << pair_info.first << pair_info.second << " "
						<< data_.size());
		PairInfos pair_infos = data_.GetEdgePairInfos(pair_info.first,
				pair_info.second);
		int cluster_index = NearestClusterIndex(pair_infos, pair_info);
		if (cluster_index >= 0) {
			MergeData(pair_infos[cluster_index], pair_info);
		} else {
			data_.AddPairInfo(pair_info, add_reversed);
		}
	}

	void RemoveEdgeInfo(EdgeId edge) {
		data_.DeleteEdgeInfo(edge);
	}

private:
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
					if (vec[j].d > min
							&& (next == -1 || vec[next].d > vec[j].d)) {
						next = j;
					}
				}
				os << vec[next].d << " " << vec[next].weight << endl;
				if (next == -1) {
					assert(false);
				}
				if (vec[next].d > 100000) {
					assert(false);
				}
				min = vec[next].d;
			}
		}
	}

	void TransferInfo(EdgeId old_edge, EdgeId new_edge, int shift = 0,
			double weight_scale = 1.0) {
		PairInfos pair_infos = GetEdgeInfo(old_edge);
		for (size_t j = 0; j < pair_infos.size(); ++j) {
			PairInfo<EdgeId> old_pair_info = pair_infos[j];
			if (old_edge != old_pair_info.second) {
				AddPairInfo(
						PairInfo<EdgeId> (new_edge, old_pair_info.second,
								old_pair_info.d - shift,
								weight_scale * old_pair_info.weight));
			} else {
				AddPairInfo(
						PairInfo<EdgeId> (new_edge, new_edge, old_pair_info.d,
								weight_scale * 0.5 * old_pair_info.weight));
			}
		}
	}

public:

//	void OutputData(ostream &os = cout) {
//		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it)
//			for (auto it1 = graph_.SmartEdgeBegin(); !it1.IsEnd(); ++it1) {
//				OutputEdgeData(*it, *it1, os);
//			}
//	}
//
//	void OutputData(string fileName) {
//		ofstream s;
//		s.open(fileName.c_str());
//		OutputData(s);
//		s.close();
//	}

	/*
	 * @return quantity of paired info
	 */
	size_t size() {
		return data_.size();
	}

	/**
	 * Method returns all data about given edge
	 */
	PairInfos GetEdgeInfo(EdgeId edge) {
		return data_.GetEdgeInfos(edge);
	}

	/**
	 * Method returns all data about distance between two given edges
	 */
	PairInfos GetEdgePairInfo(EdgeId first, EdgeId second) {
		return data_.GetEdgePairInfos(first, second);
	}

	virtual void HandleAdd(EdgeId e) {
		this->AddPairInfo(PairInfo<EdgeId> (e, e, 0, 0.0));
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}

	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
		this->AddPairInfo(PairInfo<EdgeId> (new_edge, new_edge, 0, 0.0));
		int shift = 0;
		for (size_t i = 0; i < old_edges.size(); ++i) {
			EdgeId old_edge = old_edges[i];
			TransferInfo(old_edge, new_edge, shift);
			shift -= graph_.length(old_edge);
		}
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		TransferInfo(edge1, new_edge);
		TransferInfo(edge2, new_edge);
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		double prop = (double) graph_.length(new_edge1) / graph_.length(
				old_edge);
		//		size_t shift = graph_.length(new_edge1);
		TransferInfo(old_edge, new_edge1, 0, prop);
		//		PassEdge(graph_.length(new_edge1), shift);
		TransferInfo(old_edge, new_edge2, graph_.length(new_edge1), 1 - prop);
	}

};

/**
 * This class performs the most simple offline clustering of paired information.
 */
template<class Graph>
class SimpleOfflineClusterer {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<PairInfo<EdgeId> > PairInfos;

public:
	PairedInfoIndex<Graph> &not_clustered_;

	SimpleOfflineClusterer(PairedInfoIndex<Graph> &not_clustered) :
		not_clustered_(not_clustered) {
	}

	PairInfos ProcessEdgePair(const PairInfos &infos) {
		EdgeId edge1 = infos[0].first;
		EdgeId edge2 = infos[0].second;
		double d_sum = 0;
		double weight_sum = 0;
		for (size_t i = 0; i < infos.size(); i++) {
			d_sum += infos[i].d;
			weight_sum += infos[i].weight;
		}
		PairInfo<EdgeId> sum_info(edge1, edge2, d_sum / infos.size(),
				weight_sum);
		PairInfos result;
		result.push_back(sum_info);
		return result;
	}

	void cluster(PairedInfoIndex<Graph> &clustered) {
		assert(&not_clustered_ != &clustered);
		for (typename PairedInfoIndex<Graph>::EdgePairIterator it =
				not_clustered_.begin(); it != not_clustered_.end(); ++it) {
			PairInfos newInfos = ProcessEdgePair(*it);
			for (size_t i = 0; i < newInfos.size(); i++) {
				clustered.AddPairInfo(newInfos[i]);
			}
		}
	}
};

}

#endif /* PAIRED_INFO_HPP_ */

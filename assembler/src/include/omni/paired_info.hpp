#pragma once

//#include "utils.hpp"
//#include "sequence.hpp"
//#include "paired_read.hpp"
#include "io/paired_read.hpp"
#include <cmath>
#include <map>
#include <limits>
#include <xmath.h>
#include "omni_utils.hpp"

//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace omnigraph {

/**
 * PairInfo class represents basic data unit for paired information: edges first and second appear
 * in genome at distance d_ and this information has weight weight_.
 */
template<typename EdgeId>
struct PairInfo {
	EdgeId first;
	EdgeId second;
	double d; //distance between starts. Can be negative
	double weight;
	double variance;

	PairInfo(const PairInfo& pair_info) :
			first(pair_info.first), second(pair_info.second), d(pair_info.d), weight(
					pair_info.weight), variance(pair_info.variance) {

	}

	PairInfo(EdgeId first, EdgeId second, double d, double weight,
			double variance) :
			first(first), second(second), d(d), weight(weight), variance(
					variance) {
	}

	/*	const PairInfo set_first(EdgeId first) const {
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
	 }*/

	/**
	 * Two paired infos are considered equal if they coinside in all parameters except for the weight of
	 * info.
	 */
	bool operator==(const PairInfo& rhs) const {
		const PairInfo &lhs = *this;
		return lhs.first == rhs.first && lhs.second == rhs.second && math::eq(lhs.d, rhs.d);
	}

	bool operator!=(const PairInfo& rhs) const {
		return !(*this == rhs);
	}
};

template<typename EdgeId, class Comparator>
class PairInfoComparator {
private:
	Comparator comparator_;
public:
	PairInfoComparator(Comparator comparator) : comparator_(comparator) {
	}

	bool operator()(const PairInfo<EdgeId>& lhs, const PairInfo<EdgeId>& rhs) const {
		return lhs.first == rhs.first ?
				lhs.second == rhs.second ?
						math::ls(lhs.d, rhs.d) : comparator_(lhs.second, rhs.second)
				: comparator_(lhs.first, rhs.first);
	}
};

template<typename Graph>
ostream& operator<<(ostream& os, const PairInfo<Graph>& info) {
	return os << "PairInfo: first=" << info.first << ", second=" << info.second
			<< ", distance=" << info.d << ", weight=" << info.weight
			<< ", variance=" << info.variance;
}

//typedef vector<PairInfo<> > PairInfos;

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId id) {
	return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(1)),
			-100000000/*numeric_limits<double>::min()*/, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId id) {
	return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(-1)),
			1000000000/*numeric_limits<double>::max()*/, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId e1, EdgeId e2) {
	PairInfo<EdgeId> info = MinPairInfo(e1);
	info.second = e2;
	return info;
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId e1, EdgeId e2) {
	PairInfo<EdgeId> info = MaxPairInfo(e1);
	info.second = e2;
	return info;
}

/**
 * Method returns approximate distance between occurrences of edges in genome rounded to the nearest
 * integer. In case of a tie closest to 0 value is chosen thus one can assume that distance
 * is rounded the same way as opposite one.
 * todo check that written here is true
 */
template<typename EdgeId>
int rounded_d(PairInfo<EdgeId> const& pi) {
	int res = (int) math::round(std::abs(pi.d));
	if (pi.d < 0)
		res = -res;
	return res;
}

template<typename EdgeId>
PairInfo<EdgeId> BackwardInfo(const PairInfo<EdgeId>& pi) {
	return PairInfo<EdgeId>(pi.second, pi.first, -pi.d, pi.weight, pi.variance);
}

template<typename EdgeId>
bool IsSymmetric(PairInfo<EdgeId> const& pi) {
	return pi.first == pi.second && pi.d == 0;
}

//todo try storing set<PairInfo>
template<typename EdgeId, class Comparator>
class PairInfoIndexData {
public:
	typedef set<PairInfo<EdgeId>, PairInfoComparator<EdgeId, Comparator>> Data;
	typedef typename Data::iterator data_iterator;
	typedef typename Data::const_iterator data_const_iterator;
	typedef vector<PairInfo<EdgeId>> PairInfos;

	typedef std::pair<data_const_iterator, data_const_iterator> iterator_range;


public:
	void UpdateSingleInfo(const PairInfo<EdgeId>& info, double d,
			double weight) {
		size_t count = data_.erase(info);
		VERIFY(count != 0);
		data_.insert(
				PairInfo<EdgeId>(info.first, info.second, d, weight,
						info.variance));
	}

	void ReplaceFirstEdge(const PairInfo<EdgeId>& info, EdgeId newId) {
		//		size_t count = data_.erase(info);
		//	VERIFY(count != 0);
		data_.insert(
				PairInfo<EdgeId>(newId, info.second, info.d, info.weight,
						info.variance));
	}
public:

	PairInfoIndexData(const Comparator &comparator) : comparator_(comparator), data_(PairInfoComparator<EdgeId, Comparator>(comparator_)){
	}

	data_iterator begin() const {
		auto itp = data_.begin();
		int cnt = 1;
		for (auto it = data_.begin(); it != data_.end(); ++it) {
			if (it->first == itp->first && it->second == itp->second
					&& it != data_.begin()) {
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

	data_iterator end() const {
		return data_.end();
	}

	size_t size() const {
		return data_.size();
	}

	void AddPairInfo(const PairInfo<EdgeId>& pair_info, bool addSymmetric = 1) {
//		INFO("REALLY ADD:" << pair_info.first << " " << pair_info.second << " " << pair_info.d << " " << 
//        pair_info.weight);

		data_.insert(pair_info);

		if (!IsSymmetric(pair_info) && addSymmetric)
			data_.insert(BackwardInfo(pair_info));
	}

	void DeleteEdgeInfo(EdgeId e) {
		set<PairInfo<EdgeId>, PairInfoComparator<EdgeId, Comparator>> paired_edges(comparator_);

		for (auto lower = LowerBound(e), upper = UpperBound(e); lower != upper;
				++lower) {
			paired_edges.insert(BackwardInfo(*lower));
		}

		for (auto it = paired_edges.begin(); it != paired_edges.end(); ++it) {
			data_.erase(*it);
		}

		data_.erase(LowerBound(e), UpperBound(e));
	}

	void DeletePairInfo(const PairInfo<EdgeId>& info) {

		data_.erase(info);
//		data_.erase(LowerBound(e), UpperBound(e));
	}

	void DeleteEdgePairInfo(EdgeId e1, EdgeId e2) {
		data_.erase(LowerBound(e1, e2), UpperBound(e1, e2));
		if(e1 != e2)
			data_.erase(LowerBound(e2, e1), UpperBound(e2, e1));
	}

	PairInfos GetEdgeInfos(EdgeId e) const {
		return PairInfos(LowerBound(e), UpperBound(e));
	}

	PairInfos GetEdgePairInfos(EdgeId e1, EdgeId e2) const {
		return PairInfos(LowerBound(e1, e2), UpperBound(e1, e2));
	}

	void UpdateInfo(const PairInfo<EdgeId>& info, const double d,
			const double weight, bool add_reversed) {
		UpdateSingleInfo(info, d, weight);

		if (add_reversed && !IsSymmetric(info))
			UpdateSingleInfo(BackwardInfo(info), -d, weight);
	}

	void clear() {
		data_.clear();
	}

	data_iterator LowerBound(EdgeId e) const {
		return data_.lower_bound(MinPairInfo(e));
	}

	data_iterator UpperBound(EdgeId e) const {
		return data_.upper_bound(MaxPairInfo(e));
	}

	data_iterator LowerBound(EdgeId e1, EdgeId e2) const {
		return data_.lower_bound(MinPairInfo(e1, e2));
	}

	data_iterator UpperBound(EdgeId e1, EdgeId e2) const {
		return data_.upper_bound(MaxPairInfo(e1, e2));
	}

private:
	PairInfoComparator<EdgeId, Comparator> comparator_;
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
	typedef const vector<PairInfo<EdgeId>> PairInfos;

private:
	const double max_difference_;

public:
	/**
	 * Class EdgePairIterator is used to iterate through paires of edges which have information about distance
	 * between them stored in PairedInfoIndex.
	 */
	class EdgePairIterator {
		typename PairInfoIndexData<EdgeId, typename Graph::Comparator>::data_iterator position_;
		const PairedInfoIndex<Graph> &index_;
	public:
		EdgePairIterator(
				typename PairInfoIndexData<EdgeId, typename Graph::Comparator>::data_iterator position,
				const PairedInfoIndex<Graph> &index) :
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

	EdgePairIterator begin() const {
		return EdgePairIterator(data_.begin(), *this);
	}

	EdgePairIterator end() const {
		return EdgePairIterator(data_.end(), *this);
	}

	//begin-end insert size supposed
	PairedInfoIndex(const Graph &g, double max_difference = 0.) :
			GraphActionHandler<Graph>(g, "PairedInfoIndex"), max_difference_(
					max_difference), data_(g.ReliableComparatorInstance()) {
	}

	virtual ~PairedInfoIndex() {
		TRACE("~PairedInfoIndex ok");
	}

	//	double sum() {
	//		double res = 0;
	//		for (auto it = this->g().SmartEdgeBegin(); this->g().SmartEdgeEnd() != it; ++it)
	//			for (auto it1 = this->g().SmartEdgeBegin(); this->g().SmartEdgeEnd()
	//					!= it1; ++it1) {
	//				PairInfos vec = GetEdgePairInfo(*it1, *it);
	//				if (vec.size() != 0) {
	//					for (size_t i = 0; i < vec.size(); i++)
	//						res += vec[i].weight;
	//				}
	//			}
	//		return res;
	//	}

private:

	PairInfoIndexData<EdgeId, typename Graph::Comparator> data_;

	size_t CorrectLength(const Path<EdgeId>& path, size_t index) {
		size_t result = this->g().length(path[index]);
		if (index == 0) {
			result -= path.start_pos();
		}
		if (index == path.size() - 1) {
			result -= this->g().length(path[index]) - path.end_pos();
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
		if (math::le(std::abs(info2.d - info1.d), max_difference_)) {
			return true;
		}
		return false;
	}

	void MergeData(const PairInfo<EdgeId>& info1,
			const PairInfo<EdgeId>& info2, bool add_reversed) {
		VERIFY(info1.first == info2.first && info1.second == info2.second);
		double newWeight = info1.weight + info2.weight;
		double newD = (info1.d * info1.weight + info2.d * info2.weight)
				/ newWeight;
		if (info1.first == info1.second && (info1.d == 0 || info2.d == 0)) {
			newD = 0;
		}
		data_.UpdateInfo(info1, newD, newWeight, add_reversed);
	}

	int NearestClusterIndex(const vector<PairInfo<EdgeId>>& current_pair_infos,
			const PairInfo<EdgeId>& new_info) {
		double min_dist = max_difference_;
		int answer = -1;
		for (size_t i = 0; i < current_pair_infos.size(); ++i) {
			if (math::le(std::abs(new_info.d - current_pair_infos[i].d),
					min_dist))
				if (math::ls(std::abs(new_info.d - current_pair_infos[i].d),
						min_dist) || answer == -1
						|| math::gr(std::abs(current_pair_infos[answer].d),
								std::abs(current_pair_infos[i].d))) {
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
				"IN ADD:" << pair_info.first << pair_info.second << " " << data_.size());
		PairInfos pair_infos = data_.GetEdgePairInfos(pair_info.first,
				pair_info.second);
		int cluster_index = NearestClusterIndex(pair_infos, pair_info);
		if (cluster_index >= 0) {
			MergeData(pair_infos[cluster_index], pair_info, add_reversed);
		} else {
			data_.AddPairInfo(pair_info, add_reversed);
		}
	}

	void RemoveEdgeInfo(EdgeId edge) {
		data_.DeleteEdgeInfo(edge);
	}
	void RemovePairInfo(const PairInfo<EdgeId>& pair_info) {
		data_.DeletePairInfo(pair_info);
	}

private:
	//	void OutputEdgeData(EdgeId edge1, EdgeId edge2, ostream &os = cout) {
	//		PairInfos vec = GetEdgePairInfo(edge1, edge2);
	//		if (vec.size() != 0) {
	//			os << edge1 << " " << this->g().length(edge1) << " " << edge2 << " "
	//					<< this->g().length(edge2) << endl;
	//			if (this->g().EdgeEnd(edge1) == this->g().EdgeStart(edge2))
	//				os << "+" << endl;
	//			if (this->g().EdgeEnd(edge2) == this->g().EdgeStart(edge1))
	//				os << "-" << endl;
	//			int min = INT_MIN;
	//			for (size_t i = 0; i < vec.size(); i++) {
	//				int next = -1;
	//				for (size_t j = 0; j < vec.size(); j++) {
	//					if (vec[j].d > min
	//							&& (next == -1 || vec[next].d > vec[j].d)) {
	//						next = j;
	//					}
	//				}
	//				os << vec[next].d << " " << vec[next].weight << endl;
	//				if (next == -1) {
	//					VERIFY(false);
	//				}
	//				if (vec[next].d > 100000) {
	//					VERIFY(false);
	//				}
	//				min = vec[next].d;
	//			}
	//		}
	//	}

	void TransferInfo(EdgeId old_edge, EdgeId new_edge, int shift = 0,
			double weight_scale = 1.0) {
		PairInfos pair_infos = GetEdgeInfo(old_edge);
		for (size_t j = 0; j < pair_infos.size(); ++j) {
			PairInfo<EdgeId> old_pair_info = pair_infos[j];
			if (old_edge != old_pair_info.second) {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, old_pair_info.second,
								old_pair_info.d - shift,
								weight_scale * old_pair_info.weight,
								old_pair_info.variance));
			} else if (old_pair_info.d != 0) {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, new_edge, old_pair_info.d,
								weight_scale * 0.5 * old_pair_info.weight,
								old_pair_info.variance));
			} else {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, new_edge, old_pair_info.d,
								weight_scale * old_pair_info.weight,
								old_pair_info.variance));
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
	size_t size() const {
		return data_.size();
	}

	/**
	 * Method returns all data about given edge
	 */
	PairInfos GetEdgeInfo(EdgeId edge) const {
		return data_.GetEdgeInfos(edge);
	}

	/**
	 * Method returns all data about distance between two given edges
	 */
	PairInfos GetEdgePairInfo(EdgeId first, EdgeId second) const {
		return data_.GetEdgePairInfos(first, second);
	}

	virtual void HandleAdd(EdgeId e) {
		this->AddPairInfo(PairInfo<EdgeId>(e, e, 0, 0.0, 0.));
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		this->AddPairInfo(PairInfo<EdgeId>(new_edge, new_edge, 0, 0.0, 0.));
		int shift = 0;
		for (size_t i = 0; i < old_edges.size(); ++i) {
			EdgeId old_edge = old_edges[i];
			TransferInfo(old_edge, new_edge, shift);
			shift -= this->g().length(old_edge);
		}
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		TransferInfo(edge1, new_edge);
		TransferInfo(edge2, new_edge);
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		double prop = (double) this->g().length(new_edge1)
				/ this->g().length(old_edge);
		//		size_t shift = graph_.length(new_edge1);
		TransferInfo(old_edge, new_edge1, 0, prop);
		//		PassEdge(graph_.length(new_edge1), shift);
		TransferInfo(old_edge, new_edge2, this->g().length(new_edge1),
				1 - prop);
	}

	//	bool Check() {
	//		for (auto it = data_.begin(); it != data_.end(); ++it) {
	//			PairInfo<EdgeId> inf = *it;
	//			if (inf.first != inf.second && inf.d > 0) {
	////				cout << inf.first << endl;
	////				cout << graph_.length(inf.first) << endl;
	//				if (inf.d <= graph_.length(inf.first) - 1) {
	//					cout << "gopa " << inf.first << " " << inf.second << endl;
	//					cout << graph_.length(inf.first) << " " << graph_.length(inf.second) << endl;
	//					cout << graph_.EdgeStart(inf.first) << " " << graph_.EdgeEnd(inf.first);
	//					cout << graph_.EdgeStart(inf.second) << " " << graph_.EdgeEnd(inf.second);
	//					VERIFY(false);
	//					return false;
	//				}
	//			}
	//		}
	//		return true;
	//	}

}
;

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
		VERIFY(&not_clustered_ != &clustered);
		for (typename PairedInfoIndex<Graph>::EdgePairIterator it =
				not_clustered_.begin(); it != not_clustered_.end(); ++it) {
			PairInfos newInfos = ProcessEdgePair(*it);
			for (size_t i = 0; i < newInfos.size(); i++) {
				clustered.AddPairInfo(newInfos[i]);
			}
		}
	}
};

template<size_t k, class Graph, class SequenceMapper>
class PairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef Seq<k> Kmer;

	const Graph &graph_;
	const SequenceMapper& mapper_;
	io::IReader<io::PairedRead>& stream_;
	size_t processed_count_;

	size_t CorrectLength(Path<EdgeId> path, size_t idx) {
		size_t answer = graph_.length(path[idx]);
		if (idx == 0)
			answer -= path.start_pos();
		if (idx == path.size() - 1)
			answer -= graph_.length(path[idx]) - path.end_pos();
		return answer;
	}

	void ProcessPairedRead(omnigraph::PairedInfoIndex<Graph> &paired_index,
			const io::PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();
		Path<EdgeId> path1 = mapper_.MapSequence(read1);
		Path<EdgeId> path2 = mapper_.MapSequence(read2);
		size_t distance = p_r.distance();
		int current_distance1 = distance + path1.start_pos()
				- path2.start_pos();
		for (size_t i = 0; i < path1.size(); ++i) {
			int current_distance2 = current_distance1;
			for (size_t j = 0; j < path2.size(); ++j) {
				double weight = CorrectLength(path1, i)
						* CorrectLength(path2, j);
				PairInfo<EdgeId> new_info(path1[i], path2[j], current_distance2,
						weight, 0.);
				paired_index.AddPairInfo(new_info);
				current_distance2 += graph_.length(path2[j]);
			}
			current_distance1 -= graph_.length(path1[i]);
		}
		if (++processed_count_ % 100000 == 0) {
			TRACE("Processed " << processed_count_ << " reads");
		}
	}

public:

	PairedIndexFiller(const Graph &graph, const SequenceMapper& mapper,
			io::IReader<io::PairedRead>& stream) :
			graph_(graph), mapper_(mapper), stream_(stream), processed_count_(0) {

	}

	/**
	 * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
	 */
	void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0, 0.));
		}
		INFO("Processing paired reads (takes a while)");
		stream_.reset();
		size_t n = 0;
		while (!stream_.eof()) {
			io::PairedRead p_r;
			stream_ >> p_r;
			ProcessPairedRead(paired_index, p_r);
			VERBOSE_POWER(++n, " paired reads processed");
		}
	}

private:
	DECL_LOGGER("PairedIndexFiller");

};

//New metric weight normalizer
template<class Graph>
class PairedInfoWeightNormalizer {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const size_t insert_size_;
	//todo use this param!
	const double is_var_;
	const size_t read_length_;
	const size_t k_;
	const double avg_coverage_;
public:

	//Delta better to be around 5-10% of insert size
	PairedInfoWeightNormalizer(const Graph& g, size_t insert_size, double is_var,
			size_t read_length, size_t k, double avg_coverage) :
			g_(g), insert_size_(insert_size), is_var_(is_var), read_length_(read_length),
			k_(k), avg_coverage_(avg_coverage) {
	}

	const PairInfo<EdgeId> NormalizeWeight(const PairInfo<EdgeId>& pair_info) {
		double w = 0.;
		if (math::eq(pair_info.d, 0.) && pair_info.first == pair_info.second) {
			w = 0. + g_.length(pair_info.first) - insert_size_
					+ 2 * read_length_ + 1 - k_;
		} else {
			EdgeId e1 =	(math::ge(pair_info.d, 0.)) ?
							pair_info.first : pair_info.second;
			EdgeId e2 =	(math::ge(pair_info.d, 0.)) ?
							pair_info.second : pair_info.first;
			int gap_len = std::abs(rounded_d(pair_info)) - g_.length(e1);
			int right = std::min(insert_size_,
					gap_len + g_.length(e2) + read_length_);
			int left = std::max(gap_len, int(insert_size_) - int(read_length_) - int(g_.length(e1)));
			w = 0. + right - left + 1 - k_;
		}

		double result_weight = pair_info.weight;
        if (math::gr(w, /*-10.*/0.)) {
			result_weight /= w;//(w + 10);
		}
        double cov_norm_coeff = avg_coverage_ / (2*(read_length_ - k_));
        result_weight /= cov_norm_coeff;

		PairInfo<EdgeId> result(pair_info);
		result.weight = result_weight;
		return result;
	}
};

template<class Graph>
class JumpingNormilizerFunction {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph& graph_;
	size_t read_length_;
	size_t max_norm_;

public:
	JumpingNormilizerFunction(const Graph& graph, size_t read_length, size_t max_norm) : graph_(graph), read_length_(read_length), max_norm_(max_norm) {
	}

	size_t norm(EdgeId first, EdgeId second) const {
		return std::min(std::min(graph_.length(first), graph_.length(second)), max_norm_) + read_length_ - graph_.k();
	}

	const PairInfo<EdgeId> operator()(const PairInfo<EdgeId>& pair_info) const {
		return PairInfo<EdgeId>(pair_info.first, pair_info.second, pair_info.d, pair_info.weight / norm(pair_info.first, pair_info.second), pair_info.variance);
	}
};

template<class Graph>
const PairInfo<typename Graph::EdgeId> TrivialWeightNormalization(
		const PairInfo<typename Graph::EdgeId>& pair_info) {
	return pair_info;
}

template<class Graph>
class PairedInfoNormalizer {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef boost::function<const PairInfo<EdgeId>(const PairInfo<EdgeId>&)> WeightNormalizer;
private:

	const PairedInfoIndex<Graph>& paired_index_;
	WeightNormalizer normalizing_function_;
public:

	PairedInfoNormalizer(const PairedInfoIndex<Graph>& paired_index,
			WeightNormalizer normalizing_function) :
			paired_index_(paired_index), normalizing_function_(
					normalizing_function) {

	}

	void FillNormalizedIndex(PairedInfoIndex<Graph>& normalized_index) {
		for (auto it = paired_index_.begin(); it != paired_index_.end(); ++it) {
			vector<PairInfo<EdgeId>> infos = *it;
			for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
				normalized_index.AddPairInfo(normalizing_function_(*it2));
			}
		}
	}
};

template<class Graph>
class PairedInfoSymmetryHack {
public:
	typedef typename Graph::EdgeId EdgeId;
private:
	const Graph& graph_;
	const PairedInfoIndex<Graph>& paired_index_;
public:

	PairedInfoSymmetryHack(const Graph& graph, const PairedInfoIndex<Graph>& paired_index) :
				graph_(graph), paired_index_(paired_index) {
	}

	void FillSymmetricIndex(PairedInfoIndex<Graph>& index) {
		for (auto it = paired_index_.begin(); it != paired_index_.end(); ++it) {
			vector<PairInfo<EdgeId>> infos = *it;
			for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
				index.AddPairInfo(PairInfo<EdgeId>(it2->first, it2->second, it2->d, it2->weight * 0.5, it2->variance), 0);
				index.AddPairInfo(PairInfo<EdgeId>(graph_.conjugate(it2->second), graph_.conjugate(it2->first),
						it2->d + graph_.length(it2->second) - graph_.length(it2->first), it2->weight * 0.5, it2->variance), 0);
				auto info = index.GetEdgePairInfo(infos[0].first, infos[0].second);
				auto symmetric_info = index.GetEdgePairInfo(graph_.conjugate(info[0].second), graph_.conjugate(info[0].first));
			}
		}
	}
};


}

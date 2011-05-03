#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "strobe_reader.hpp"
#include "sequence.hpp"
#include <map>

#define MERGE_DATA_ABSOLUTE_DIFFERENCE 1000
//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace de_bruijn {

template<class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	template<size_t kmer_size>
	void ProcessPairedRead(PairedRead &p_r,
			de_bruijn::SimpleSequenceMapper<kmer_size, Graph> &read_threader) {
		Sequence read1 = p_r.first().getSequence();
		Sequence read2 = p_r.second().getSequence();
		de_bruijn::Path<EdgeId> path1 = read_threader.MapSequence(read1);
		de_bruijn::Path<EdgeId> path2 = read_threader.MapSequence(read2);
		size_t distance = CountDistance(read1, read2);
		int current_distance1 = distance + path1.start_pos()
				- path2.start_pos();
		for (size_t i = 0; i < path1.size(); ++i) {
			int current_distance2 = current_distance1;
			for (size_t j = 0; j < path2.size(); ++j) {
				double weight = CorrectLength(path1, i) * CorrectLength(path2,
						j);
				PairInfo
						new_info(path1[i], path2[j], current_distance2, weight);
				AddPairInfo(new_info);
				current_distance2 += graph_.length(path2[j]);
			}
			current_distance1 -= graph_.length(path1[i]);
		}
	}
public:

	//begin-end insert size supposed
	PairedInfoIndex(Graph &g, size_t insert_size) :
		graph_(g), insert_size_(insert_size) {
		g.AddActionHandler(this);
	}

	virtual ~PairedInfoIndex() {
		graph_.RemoveActionHandler(this);
	}

	template<size_t kmer_size, class Stream>
	void FillIndex(const DeBruijnPlus<kmer_size + 1, EdgeId>& index,
			Stream& stream) {
		data_.clear();
		de_bruijn::SmartEdgeIterator<Graph> it = graph_.SmartEdgeBegin();
		for (de_bruijn::SmartEdgeIterator<Graph> it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it; ++it) {
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
		for (de_bruijn::SmartEdgeIterator<Graph> it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd()
				!= it; ++it)
			for (de_bruijn::SmartEdgeIterator<Graph> it1 =
					graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it1; ++it1) {
				PairInfos vec = GetEdgePairInfo(*it1, *it);
				if (vec.size() != 0) {
					for (size_t i = 0; i < vec.size(); i++)
						res += vec[i].weight_;
				}
			}
		return res;
	}

	class PairInfo {
		friend class PairedInfoIndex;
		friend class PairedInfoIndexData;
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
		typedef multimap<pair<EdgeId, EdgeId> , pair<double, double>> Data;
		typedef typename Data::iterator data_iterator;
		typedef typename Data::const_iterator const_data_iterator;

		Data data_;

		data_iterator LowerBound(EdgeId e) {
			return data_.lower_bound(make_pair(e, (EdgeId) 0));
		}

		data_iterator UpperBound(EdgeId e) {
			return data_.upper_bound(make_pair(e, (EdgeId) ((size_t) -1)));
		}

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
	};

	Graph& graph_;
	PairInfoIndexData data_;

	//begin-end insert size supposed
	size_t insert_size_;

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

	size_t CountDistance(const Sequence& read1, const Sequence& read2) {
		return insert_size_ - read2.size();
	}

	bool CanMergeData(const PairInfo& info1, const PairInfo& info2) {
		if (info1.first_ != info2.first_ || info1.second_ != info2.second_)
			return false;
		if (std::abs(info2.d_ - info1.d_) <= MERGE_DATA_ABSOLUTE_DIFFERENCE) {
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

	void RemoveEdgeInfo(EdgeId edge) {
		data_.DeleteEdgeInfo(edge);
	}

	void OutputEdgeData(EdgeId edge1, EdgeId edge2) {
		PairInfos vec = GetEdgePairInfo(edge1, edge2);
		if (vec.size() != 0) {
			cout << edge1 << " " << graph_.length(edge1) << " " << edge2 << " "
					<< graph_.length(edge2) << endl;
			if (graph_.EdgeEnd(edge1) == graph_.EdgeStart(edge2))
				cout << "+" << endl;
			if (graph_.EdgeEnd(edge2) == graph_.EdgeStart(edge1))
				cout << "-" << endl;
			int min = INT_MIN;
			for (size_t i = 0; i < vec.size(); i++) {
				int next = -1;
				for (size_t j = 0; j < vec.size(); j++) {
					if (vec[j].d() > min && (next == -1 || vec[next].d()
							> vec[j].d())) {
						next = j;
					}
				}
				cout << vec[next].d() << " " << vec[next].weight() << endl;
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

	void OutputData() {
		for (de_bruijn::SmartEdgeIterator<Graph> it = graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd()
				!= it; ++it)
			for (de_bruijn::SmartEdgeIterator<Graph> it1 =
					graph_.SmartEdgeBegin(); graph_.SmartEdgeEnd() != it1; ++it1) {
				OutputEdgeData(*it, *it1);
			}
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
}

#endif /* PAIRED_INFO_HPP_ */

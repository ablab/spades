#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "strobe_reader.hpp"
#include "sequence.hpp"
#include <map>

#define MERGE_DATA_ABSOLUTE_DIFFERENCE 3
#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace de_bruijn {

template<size_t kmer_size, class Stream, class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	struct PairInfo {
		EdgeId first_;
		EdgeId second_;
		int d_;//distance between starts. Can be negative
		double weight_;
		PairInfo(EdgeId first, EdgeId second, int d, double weight) :
			first_(first), second_(second), d_(d), weight_(weight) {
		}

		bool operator==(const PairInfo& rhs) const {
			return first_ == rhs.first_ && second_ == rhs.second_ && d_
					== rhs.d_/* && weight_ == rhs.weight_*/;
		}
	};

	typedef const vector<PairInfo> PairInfos;
	typedef typename PairInfos::const_iterator infos_iterator;

	//	template<size_t kmer_size>
	PairedInfoIndex(Graph &g, const SimpleIndex<kmer_size + 1, EdgeId>& index,
			Stream stream) :
		graph_(g) {
		CollectData/*<kmer_size> */(index, stream);
		g.AddActionHandler(this);
	}

private:
	//todo try storing set<PairInfo>
	class PairInfoIndexData {
		typedef multimap<pair<EdgeId, EdgeId> , pair<size_t, double>> Data;
		typedef typename Data::iterator data_iterator;
		typedef typename Data::const_iterator const_data_iterator;

		Data data_;

		data_iterator LowerBound(EdgeId e) {
			return data_.lower_bound(make_pair(e, (EdgeId)0));
		}

		data_iterator UpperBound(EdgeId e) {
			return data_.upper_bound(make_pair(e, (EdgeId)((size_t) -1)));
		}

		const PairInfo BackwardInfo(const PairInfo& pair_info) {
			return PairInfo(pair_info.second_, pair_info.first_, -pair_info.d_,
					pair_info.weight_);
		}

		const PairInfo AsPairInfo(
				const pair<pair<EdgeId, EdgeId> , pair<size_t, double>>& pair) {
			return PairInfo(pair.first.first, pair.first.second,
					pair.second.first, pair.second.second);
		}

		const pair<EdgeId, EdgeId> EdgePair(const PairInfo& pair_info) {
			return make_pair(pair_info.first_, pair_info.second_);
		}

		const pair<pair<EdgeId, EdgeId> , pair<size_t, double>> AsPairOfPairs(
				const PairInfo& pair_info) {
			return make_pair(EdgePair(pair_info),
					make_pair(pair_info.d_, pair_info.weight_));
		}

		void UpdateSingleInfo(const PairInfo& info, const int d,
				const double weight) {
			bool updated = false;
			for (typename Data::iterator lower = data_.lower_bound(EdgePair(info)),
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
			data_.insert(AsPairOfPairs(pair_info));
			data_.insert(AsPairOfPairs(BackwardInfo(pair_info)));
		}

		void DeleteEdgeInfo(EdgeId e) {
			set<EdgeId> paired_edges;
			for (const_data_iterator lower = LowerBound(e), upper = UpperBound(e); lower
					!= upper; ++lower) {
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
			for (const_data_iterator lower = LowerBound(e), upper = UpperBound(e); lower
					!= upper; ++lower) {
				answer.push_back(AsPairInfo(*lower));
			}
			return answer;
		}

		PairInfos GetEdgePairInfos(EdgeId e1, EdgeId e2) {
			vector<PairInfo> answer;
			for (const_data_iterator lower = data_.lower_bound(make_pair(e1, e2)),
					upper = data_.upper_bound(make_pair(e1, e2)); lower != upper; ++lower) {
				answer.push_back(AsPairInfo(*lower));
			}
			return answer;
		}

		void UpdateInfo(const PairInfo& info, const int d, const double weight) {
			UpdateSingleInfo(info, d, weight);
			UpdateSingleInfo(BackwardInfo(info), -d, weight);
		}

	};

	Graph& graph_;
	PairInfoIndexData data_;

	size_t CorrectLength(const de_bruijn::Path<EdgeId>& path, size_t index) {
		if (index == 0) {
			return graph_.length(path[index]) - path.start_pos();
		}
		if (index == path.size() - 1) {
			return path.end_pos();
		}
		return graph_.length(path[index]);
	}

	void PassEdge(size_t edge_length, size_t &path_nucls_passed) {
		if (path_nucls_passed == 0) {
			path_nucls_passed += kmer_size;
		}
		path_nucls_passed += edge_length;
	}

	//	template<size_t kmer_size>
	void CollectData(const SimpleIndex<kmer_size + 1, EdgeId>& index,
			Stream &stream) {
		//todo
		size_t d = 100;

		typedef Seq<kmer_size + 1> KPOMer;
		de_bruijn::SimpleReadThreader<kmer_size, Graph> read_threader(graph_,
				index);
		while (!stream.eof()) {
			vector<Read> reads;
			stream >> reads;
			de_bruijn::Path<EdgeId> path1 = read_threader.ThreadRead(
					reads[0].getSequence());
			de_bruijn::Path<EdgeId> path2 = read_threader.ThreadRead(
					reads[1].getSequence());
			//walken path lengths
			size_t length1 = 0;
			size_t length2 = 0;
			for (size_t i = 0; i < path1.size(); ++i) {
				for (size_t j = 0; j < path2.size(); ++j) {
					AddPairInfo(PairInfo(path1[i], path2[j], d + length2 - length1,
							CorrectLength(path1, i) * CorrectLength(path2, j)));
					PassEdge(CorrectLength(path2, j), length2);
				}
				PassEdge(CorrectLength(path1, i), length1);
			}
		}
	}

	bool MergeData(const PairInfo& info1, const PairInfo& info2) {
		if (info1.first_ != info2.first_ || info1.second_ != info2.second_)
			return false;

		if (std::abs(info2.d_ - info1.d_) <= MERGE_DATA_ABSOLUTE_DIFFERENCE
				&& std::abs(info2.d_ - info1.d_) <= info1.d_
						* MERGE_DATA_RELATIVE_DIFFERENCE) {
			double newWeight = info1.weight_ + info2.weight_;
			int newD = std::floor(
					(info1.d_ * info1.weight_ + info2.d_ * info2.weight_) / info2.weight_ + 0.5);
			data_.UpdateInfo(info1, newD, newWeight);
			return true;
		}
		return false;
	}

//	bool MergeData(PairInfo info, const int d, const double weight) {
//		if (std::abs(d - info.d_) <= MERGE_DATA_ABSOLUTE_DIFFERENCE
//				&& std::abs(d - info.d_) <= info.d_
//						* MERGE_DATA_RELATIVE_DIFFERENCE) {
//			double newWeight = info.weight_ + weight;
//			int newD = std::floor(
//					(info.d_ * info.weight_ + d * weight) / weight + 0.5);
//			data_.UpdateInfo(info, newD, newWeight);
//			return true;
//		}
//		return false;
//	}

	void AddPairInfo(const PairInfo& pair_info/*const EdgeId first, const EdgeId second, const int d,
			const double weight*/) {
		PairInfos pair_infos = data_.GetEdgePairInfos(pair_info.first_, pair_info.second_);
		for (infos_iterator it = pair_infos.begin(); it != pair_infos.end(); ++it) {
			if (MergeData(*it, pair_info)) {
				return;
			}
		}
		data_.AddPairInfo(pair_info);
	}

	void RemoveEdgeInfo(EdgeId edge) {
		data_.DeleteEdgeInfo(edge);
	}

	void TransferInfo(EdgeId old_edge, EdgeId new_edge, size_t shift = 0, double weight_scale = 1.0) {
		PairInfos pair_infos = GetEdgeInfo(old_edge);
		for (size_t j = 0; j < pair_infos.size(); ++j) {
			PairInfo old_pair_info = pair_infos[j];
			AddPairInfo(PairInfo(new_edge, old_pair_info.second_, old_pair_info.d_ - shift, weight_scale * old_pair_info.weight_));
		}
		RemoveEdgeInfo(old_edge);
	}

public:
	PairInfos GetEdgeInfo(EdgeId edge) {
		return data_.GetEdgeInfos(edge);
	}

	PairInfos GetEdgePairInfo(EdgeId first, EdgeId second) {
		return data_.GetEdgePairInfos(first, second);
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}


	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
		size_t shift = 0;
		for (size_t i = 0; i < old_edges.size(); ++i) {
			EdgeId old_edge = old_edges[i];
			TransferInfo(old_edge, new_edge, shift);
			PassEdge(graph_.length(old_edge), shift);
		}
	}

	virtual void HandleGlue(EdgeId old_edge, EdgeId new_edge) {
		TransferInfo(old_edge, new_edge);
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		double prop = (double) graph_.length(new_edge1) / graph_.length(old_edge);
		size_t shift = 0;
		TransferInfo(old_edge, new_edge1, shift, prop);
		PassEdge(graph_.length(new_edge1), shift);
		TransferInfo(old_edge, new_edge1, shift, 1-prop);
	}

	virtual ~PairedInfoIndex() {
		graph_.RemoveActionHandler(this);
	}

};
}

#endif /* PAIRED_INFO_HPP_ */

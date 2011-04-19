#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "strobe_reader.hpp"
#include "sequence.hpp"

#define MERGE_DATA_ABSOLUTE_DIFFERENCE 3
#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace de_bruijn {

template<class Graph>
class PairedInfoIndex: GraphActionHandler<Graph> {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	struct PairInfo {
		const EdgeId first_;
		const EdgeId second_;
		int d_;//distance between starts. Can be negative
		double weight_;
		PairInfo(EdgeId first, EdgeId second, int d, double weight) :
			first_(first), second_(second), d_(d), weight_(weight) {
		}
	};

private:
	Graph &graph_;
	map<EdgeId, vector<PairInfo> > data_;

	template<size_t kmer_size>
	void collectData(SimpleIndex<kmer_size, EdgeId> &index,
			StrobeReader<2, Read, ireadstream> &reader) {
		//TODO
	}

	void UpdateInfo(PairInfo &info, const int d, const double weight) {
		for (typename vector<PairInfo>::iterator it =
				data_[info.second_].begin(); it != data_[info.second_].end(); ++it) {
			if (it->second_ == info.first_ && it->d_ == info.d_ && it->weight_
					== info.weight_) {
				it->d_ = info.d_;
				it->weight_ = info.weight_;
				break;
			}
		}
		info.d_ = d;
		info.weight_ = weight;
	}

	bool MergeData(PairInfo info, const int d, const double weight) {
		if (std::abs(d - info.d_) <= MERGE_DATA_ABSOLUTE_DIFFERENCE
				&& std::abs(d - info.d_) <= info.d_
						* MERGE_DATA_RELATIVE_DIFFERENCE) {
			double newWeight_ = info.weight_ + weight;
			int newD = std::floor(
					(info.d_ * info.weight + d * weight) / weight + 0.5);
			update(info, d, weight);
			return true;
		}
		return false;
	}

	void AddPairInfoToData(const EdgeId first, const EdgeId second,
			const int d, const double weight) {
		data_[first].push_back(PairData(first, second, d, weight));
		data_[second].push_back(PairData(second, first, -d, weight));
	}

	void AddPairInfo(const EdgeId first, const EdgeId second, const int d,
			const double weight) {
		typename map<EdgeId, vector<PairInfo> >::iterator it =
				data_.find(first);
		if (it == data_.end()) {
			AddPairInfoToData(first, second, d, weight);
		} else {
			vector<PairInfo> &edgeData = *it;
			for (size_t i = 0; i < edgeData.size(); i++)
				if (edgeData[i].second == second) {
					if (MergeData(edgeData[i], d, weight)) {
						break;
					}
				}
		}
	}

	void RemoveEdgeInfo(const EdgeId edge) {
		data_.erase(edge);
	}

public:
	template<size_t kmer_size>
	PairedInfoIndex(Graph &g, SimpleIndex<kmer_size, EdgeId> &index,
			StrobeReader<2, Read, ireadstream> &reader) :
		graph_(g) {
		collectData<kmer_size> (index, reader);
		g.AddActionHandler(this);
	}

	vector<const PairInfo> GetEdgeInfo(EdgeId edge) {
		typename map<EdgeId, vector<PairInfo> >::iterator it = data_.find(edge);
		if (it == data_.end()) {
			vector<const PairInfo> res;
			return res;
		} else {
			return vector<const PairInfo> (it->begin(), it->end());
		}
	}

	vector<const PairInfo> GetEdgePairInfo(EdgeId first, EdgeId second) {
		typename map<EdgeId, vector<PairInfo> >::iterator it =
				data_.find(first);
		vector<const PairInfo> res;
		if (it == data_.end()) {
			return res;
		} else {
			for (typename vector<PairInfo>::iterator dataIt = it->begin(); it
					!= it->end(); ++it) {
				assert(dataIt->first == first);
				if (dataIt->second_ == second) {
					res.push_back(*it);
				}
			}
			return res;
		}
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}

	virtual void HandleMerge(vector<EdgeId> oldEdge, EdgeId newEdge) {
		//TODO
	}

	virtual void HandleGlue(EdgeId oldEdge, EdgeId newEdge) {
		//TODO
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		//TODO
	}

	~PairedInfoIndex() {
		graph_.RemoveActionHandler(this);
	}

};
}

#endif /* PAIRED_INFO_HPP_ */

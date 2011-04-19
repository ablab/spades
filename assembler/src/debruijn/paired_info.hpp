#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "strobe_reader.hpp"
#include "sequence.hpp"
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

	void AddPairInfo(const EdgeId first, const EdgeId second, const int d_,
			const double weight) {
		//TODO
	}

	void RemoveEdgeInfo(const EdgeId edge) {
		//TODO
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
		//TODO
	}

	vector<const PairInfo> GetEdgePairInfo(EdgeId first, EdgeId second) {
		//TODO
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

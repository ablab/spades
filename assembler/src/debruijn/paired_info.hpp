#ifndef PAIRED_INFO_HPP_
#define PAIRED_INFO_HPP_
#include "utils.hpp"
#include "strobe_reader.hpp"
#include "sequence.hpp"
namespace de_bruijn {

template<size_t kmer_size, class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {
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

//	template<size_t kmer_size>
	PairedInfoIndex(Graph &g, const SimpleIndex<kmer_size + 1, EdgeId>& index,
			StrobeReader<2, Read, ireadstream> &reader) :
		graph_(g) {
		CollectData/*<kmer_size> */(index, reader);
		g.AddActionHandler(this);
	}

private:
	Graph& graph_;
	map<EdgeId, vector<PairInfo> > data_;

	size_t CorrectLength(const de_bruijn::Path<EdgeId>& path, size_t index) {
		if (index == 0) {
			return graph_.length(path[index]) - path.start_pos();
		}
		if (index == path.size() - 1) {
			return path.end_pos();
		}
		return graph_.length(path[index]);
	}

//	template<size_t kmer_size>
	void CollectData(const SimpleIndex<kmer_size + 1, EdgeId>& index,
			StrobeReader<2, Read, ireadstream> &reader) {
		//todo
		size_t d = 100;

		typedef Seq<kmer_size + 1> KPOMer;
		de_bruijn::SimpleReadThreader<kmer_size, Graph> read_threader(graph_, index);
		while (!reader.eof()) {
			vector<Read> reads;
			reader >> reads;
			de_bruijn::Path<EdgeId> path1 = read_threader.ThreadRead(reads[0].getSequence());
			de_bruijn::Path<EdgeId> path2 = read_threader.ThreadRead(reads[1].getSequence());
			//walken path lengths
			size_t length1 = 0;
			size_t length2 = 0;
			for (size_t i = 0; i < path1.size(); ++i) {
				for (size_t j = 0; j < path2.size(); ++j) {
					AddPairInfo(path1[i], path2[j], d + length2 - length1, CorrectLength(path1, i) * CorrectLength(path2, j));
					if (length2 == 0) {
						length2 += kmer_size;
					}
					length2 += CorrectLength(path2, j);
				}
				if(length1 == 0) {
					length1 += kmer_size;
				}
				length1 += CorrectLength(path1, i);
			}
		}
	}

	void AddPairInfo(const EdgeId first, const EdgeId second, const int d_,
			const double weight) {
		//TODO
	}

	void RemoveEdgeInfo(const EdgeId edge) {
		//TODO
	}

public:
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

	virtual ~PairedInfoIndex() {
		graph_.RemoveActionHandler(this);
	}

};
}

#endif /* PAIRED_INFO_HPP_ */

#ifndef OMNI_TOOLS_HPP_
#define OMNI_TOOLS_HPP_

#include "omni_utils.hpp"
#include "paired_info.hpp"

namespace omnigraph {

/**
 * Compresser compresses vertices with unique incoming and unique outgoing edge in linear time while
 * simple one-by-one compressing has square complexity.
 */
template<class Graph>
class Compressor {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;

	bool GoUniqueWayForward(EdgeId &e) {
		VertexId u = graph_.EdgeEnd(e);
		if (!graph_.CheckUniqueOutgoingEdge(u)
				|| !graph_.CheckUniqueIncomingEdge(u)) {
			return false;
		}
		e = graph_.GetUniqueOutgoingEdge(u);
		return true;
	}

	bool GoUniqueWayBackward(EdgeId &e) {
		VertexId u = graph_.EdgeStart(e);
		if (!graph_.CheckUniqueOutgoingEdge(u)
				|| !graph_.CheckUniqueIncomingEdge(u)) {
			return false;
		}
		e = graph_.GetUniqueIncomingEdge(u);
		return true;
	}

public:
	Compressor(Graph &graph) :
		graph_(graph) {
	}

	/**
	 * Method compresses longest possible path, containing given vertex.
	 * @param vertex to be compressed as part of a path
	 * @return true if vertex can be compressed and false otherwise
	 */
	bool CompressVertex(VertexId v) {
		TRACE("Processing vertex " << v << " started");
		if (!graph_.CheckUniqueOutgoingEdge(v)
				|| !graph_.CheckUniqueIncomingEdge(v)) {
			TRACE("Vertex " << v << " judged NOT compressible. Proceeding to the next vertex");
			TRACE("Processing vertex " << v << " finished");
			return false;
		}
		TRACE("Vertex " << v << " judged compressible");
		EdgeId e = graph_.GetUniqueOutgoingEdge(v);
		EdgeId start_edge = e;
		while (GoUniqueWayBackward(e) && e != start_edge) {
		}
		vector<EdgeId> mergeList;
		//		e = graph_.conjugate(e);
		start_edge = e;
		do {
			mergeList.push_back(e);
		} while (GoUniqueWayForward(e) && e != start_edge);
		EdgeId new_edge = graph_.MergePath(mergeList);
		TRACE("Vertex " << v << " compressed and is now part of edge " << new_edge);
		TRACE("Processing vertex " << v << " finished");
		return true;
	}

	/**
	 * Method compresses all vertices which can be compressed.
	 */
	void CompressAllVertices() {
		TRACE("Vertex compressing started");
		//SmartVertexIterator<Graph> end = graph_.SmartVertexEnd();
		for (auto it = graph_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			VertexId v = *it;
			CompressVertex(v);
		}
		TRACE("Vertex compressing finished")
	}

private:
	DECL_LOGGER("Compressor")
};

template<class Graph>
class Cleaner {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;

public:
	Cleaner(Graph &graph) :
		graph_(graph) {
	}

	void Clean() {
		for (auto iter = graph_.SmartVertexBegin(); !iter.IsEnd(); ++iter) {
			if (graph_.IsDeadStart(*iter) && graph_.IsDeadEnd(*iter)) {
				graph_.DeleteVertex(*iter);
			}
		}
	}

private:
	DECL_LOGGER("Compressor")
};

template<class Graph>
class GraphCopier {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph &graph_;
public:
	GraphCopier(Graph &graph) :
		graph_(graph) {
	}
	template<class CopyGraph>
	void StructureCopy(CopyGraph &new_graph,
			map<VertexId, typename CopyGraph::VertexId> &VerticesCopies,
			map<EdgeId, typename CopyGraph::EdgeId> &EdgesCopies) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
			NewVertexId new_vertex = new_graph.AddVertex(graph_.data(*iter));
			DEBUG("Added vertex "<< new_vertex);
			VerticesCopies.insert(make_pair(*iter, new_vertex));
		}
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			EdgeId edge = *iter;
			NewEdgeId new_edge = new_graph.AddEdge(
					VerticesCopies[graph_.EdgeStart(edge)],
					VerticesCopies[graph_.EdgeEnd(edge)], graph_.data(edge));
			EdgesCopies.insert(make_pair(*iter, new_edge));
		}
	}
	template<class CopyGraph>
	void Copy(CopyGraph &new_graph) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		map<EdgeId, NewEdgeId> EdgesCopies;
		map<VertexId, NewVertexId> VerticesCopies;
		StructureCopy<CopyGraph> (new_graph, VerticesCopies, EdgesCopies);
	}
	template<class CopyGraph>
	void TotalCopy(CopyGraph &new_graph) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		map<EdgeId, NewEdgeId> EdgesCopies;
		map<VertexId, NewVertexId> VerticesCopies;
		StructureCopy<CopyGraph> (new_graph, VerticesCopies, EdgesCopies);
	}
};

template<class Graph>
class ConjugateGraphCopier {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph &graph_;
public:
	ConjugateGraphCopier(Graph &graph) :
		graph_(graph) {
	}

	template<class CopyGraph>
	void Copy(CopyGraph new_graph) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		map<VertexId, NewVertexId> copy;
		for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
			if (copy.count(*iter) == 0) {
				NewVertexId new_vertex =
						new_graph.AddVertex(graph_.data(*iter));
				copy.insert(make_pair(*iter, new_vertex));
				copy.insert(
						make_pair(graph_.conjugate(*iter),
								new_graph.conjugate(new_vertex)));
			}
		}
		set<EdgeId> was;
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (was.count(*iter) == 0) {
				new_graph.AddEdge(copy[graph_.EdgeStart(*iter)],
						copy[graph_.EdgeEnd(*iter)], graph_.data(*iter));
				was.insert(*iter);
				was.insert(graph_.conjugate(*iter));
			}
		}
	}
};

template<class Graph>
class TrivialEdgePairChecker {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	const size_t bound_;
public:
	TrivialEdgePairChecker(Graph &graph, size_t bound = (size_t) -1) :
		graph_(graph), bound_(bound) {
	}

	/*
	 * Very bad code. Shame on me.
	 */
	bool GoForward(EdgeId &edge) {
		if (!graph_.CheckUniqueOutgoingEdge(graph_.EdgeEnd(edge))) {
			return false;
		}
		edge = graph_.GetUniqueOutgoingEdge(graph_.EdgeEnd(edge));
		return true;
	}

	bool GoBackward(EdgeId &edge) {
		if (!graph_.CheckUniqueIncomingEdge(graph_.EdgeStart(edge))) {
			return false;
		}
		edge = graph_.GetUniqueIncomingEdge(graph_.EdgeStart(edge));
		return true;
	}

	bool CheckForward(EdgeId edge1, EdgeId edge2) {
		set<EdgeId> was;
		size_t length = 0;
		do {
			if (edge1 == edge2)
				return true;
			if (was.count(edge1) != 0)
				return false;
			was.insert(edge1);
			length += graph_.length(edge1);
		} while (length <= bound_ && GoForward(edge1));
		return false;
	}

	bool CheckBackward(EdgeId edge1, EdgeId edge2) {
		set<EdgeId> was;
		size_t length = 0;
		do {
			if (edge1 == edge2)
				return true;
			if (was.count(edge1) != 0)
				return false;
			was.insert(edge1);
			length += graph_.length(edge1);
		} while (length <= bound_ && GoBackward(edge1));
		return false;
	}

	bool Check(EdgeId edge1, EdgeId edge2) {
		return CheckForward(edge1, edge2) || CheckBackward(edge2, edge1)
		/*|| CheckForward(edge2, edge1) || CheckBackward(edge1, edge2)*/;
	}
};

template<class Graph>
class PairInfoFilter {
private:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	double weight_threshold_;

	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
		double s = 0.0;
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if ((*it).d > graph_.length((*it).first)) {
				s += (*it).weight;
			}
		}
		return s > weight_threshold_;
	}

	void InsertData(const vector<PairInfo<EdgeId>>& infos,
			PairedInfoIndex<Graph> &new_index) {
		double s = 0.0;
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if ((*it).d > graph_.length((*it).first)) {
				s += (*it).weight;
				new_index.AddPairInfo(*it);
			}
		}
	}

public:
	PairInfoFilter(Graph &graph, double weight_threshold) :
		graph_(graph), weight_threshold_(weight_threshold) {

	}

	void Filter(PairedInfoIndex<Graph> &old_index,
			PairedInfoIndex<Graph> &new_index) {
		for (auto iterator = old_index.begin(); iterator != old_index.end(); ++iterator) {
			auto data = *iterator;
			InsertData(data, new_index);
		}
	}
};

}

#endif /* OMNI_TOOLS_HPP_ */

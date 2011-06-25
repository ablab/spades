#ifndef OMNI_TOOLS_HPP_
#define OMNI_TOOLS_HPP_

#include "omni_utils.hpp"

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
	void Copy(CopyGraph new_graph) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		map<VertexId, NewVertexId> copy;
		for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
			NewVertexId new_vertex = new_graph.AddVertex(graph_.data(*iter));
			copy.insert(make_pair(*iter, new_vertex));
		}
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			EdgeId edge = *iter;
			new_graph.AddEdge(copy[graph_.EdgeStart(edge)],
					copy[graph_.EdgeEnd(edge)], graph_.data(edge));
		}
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
				NewVertexId new_vertex = new_graph.AddVertex(
						graph_.data(*iter));
				copy.insert(make_pair(*iter, new_vertex));
				copy.insert(
						make_pair(graph_.conjugate(*iter),
								new_graph.conjugate(new_vertex)));
			}
		}
		set<EdgeId> was;
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (was.count(*iter) == 0) {
				new_graph.AddEdge(copy[graph_.EdgeStart(*iter)], copy[graph_.EdgeEnd(*iter)],
						graph_.data(*iter));
				was.insert(*iter);
				was.insert(graph_.conjugate(*iter));
			}
		}
	}
};

}

#endif /* OMNI_TOOLS_HPP_ */

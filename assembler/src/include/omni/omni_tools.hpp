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
		while (GoUniqueWayBackward(e) && e != start_edge
				&& !graph_.RelatedVertices(graph_.EdgeStart(e), graph_.EdgeEnd(e))) {
		}
		vector<EdgeId> mergeList;
		//		e = graph_.conjugate(e);
		start_edge = e;
		do {
			mergeList.push_back(e);
		} while (GoUniqueWayForward(e) && e != start_edge
				&& !graph_.RelatedVertices(graph_.EdgeStart(e), graph_.EdgeEnd(e)));
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

		//in current implementation will work incorrectly if smart iterator won't give vertex and its conjugate
		//(in case of self-conjugate edges)
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
class IsolatedEdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	size_t max_length_;

	bool IsTerminalVertex(VertexId v) {
		return g_.IncomingEdgeCount(v) + g_.OutgoingEdgeCount(v) == 1;
	}

public:
	IsolatedEdgeRemover(Graph& g, size_t max_length): g_(g), max_length_(max_length) {
	}

	void RemoveIsolatedEdges() {
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (IsTerminalVertex(g_.EdgeStart(*it)) && IsTerminalVertex(g_.EdgeEnd(*it))) {
				g_.DeleteEdge(*it);
			}
		}
		Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
	}

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
	void Copy(CopyGraph& new_graph) {
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
	const Graph &graph_;
	const size_t bound_;
public:
	TrivialEdgePairChecker(const Graph &graph, size_t bound = (size_t) -1) :
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
	const Graph &graph_;
	double weight_threshold_;

//	bool ContainsPositiveDistance(const vector<PairInfo<EdgeId>>& infos) {
//		double s = 0.0;
//		for (auto it = infos.begin(); it != infos.end(); ++it) {
//			if ((*it).d > graph_.length((*it).first)) {
//				s += (*it).weight;
//			}
//		}
//		return s > weight_threshold_;
//	}
//
//	void InsertData(const vector<PairInfo<EdgeId>>& infos,
//			PairedInfoIndex<Graph> &new_index) {
//		double s = 0.0;
//		for (auto it = infos.begin(); it != infos.end(); ++it) {
//			if ((*it).d > graph_.length((*it).first)) {
//				s += (*it).weight;
//				new_index.AddPairInfo(*it);
//			}
//		}
//	}

public:
	PairInfoFilter(const Graph &graph, double weight_threshold) :
		graph_(graph), weight_threshold_(weight_threshold) {

	}

	void Filter(const PairedInfoIndex<Graph> &old_index,
			PairedInfoIndex<Graph> &new_index) {
		for (auto it = old_index.begin(); it != old_index.end(); ++it) {
			auto infos = *it;
			for (auto info_it = infos.begin(); info_it != infos.end(); ++info_it) {
				if (math::ge(info_it->weight, weight_threshold_)) {
					new_index.AddPairInfo(*info_it, false);
				}
			}
		}
	}
};

template<class Graph>
class PairInfoChecker {
private:
	typedef typename Graph::EdgeId EdgeId;
	const EdgesPositionHandler<Graph> &positions_;
	size_t first_bound_;
	const size_t second_bound_;
	vector<double> perfect_matches_;
	vector<double> good_matches_;
	vector<double> mismatches_;
	vector<double> imperfect_matches_;

public:
	PairInfoChecker(const EdgesPositionHandler<Graph> &positions,
		size_t first_bound, size_t second_bound) :
			positions_(positions), first_bound_(first_bound), second_bound_(
					second_bound) {
	}

	void Check(const PairedInfoIndex<Graph> &paired_index) {
		for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
			auto vec = *it;
			for (auto vec_it = vec.begin(); vec_it != vec.end(); ++vec_it) {
				size_t code = CheckSingleInfo(*vec_it);
				if (code == 0) {
					perfect_matches_.push_back(vec_it->weight);
				} else if (code == 1) {
					good_matches_.push_back(vec_it->weight);
				} else if (code == 2) {
					mismatches_.push_back(vec_it->weight);
				} else if (code == 3) {
					imperfect_matches_.push_back(vec_it->weight);
				}
			}
		}
	}

	size_t CheckSingleInfo(PairInfo<EdgeId> info) {
		const vector<EdgePosition> &pos1 = positions_.GetEdgePositions(info.first);
		const vector<EdgePosition> &pos2 = positions_.GetEdgePositions(info.second);
		bool good_match_found = false;
		for (size_t i = 0; i < pos1.size(); i++)
			for (size_t j = 0; j < pos2.size(); j++) {
				if (abs(pos1[i].start_ + info.d - pos2[j].start_)
						<= first_bound_ + info.variance) {
					if (info.variance == 0) {
						return 0;
					} else {
						return 3;
					}
				} else if (abs(pos1[i].start_ + info.d - pos2[j].start_)
						<= second_bound_) {
					good_match_found = true;
				}
			}
		if(good_match_found) {
			return 1;
		} else {
			return 2;
		}
	}

	void WriteResultsToFile(vector<double> results, const string &file_name) {
		sort(results.begin(), results.end());
		ofstream os;
		os.open(file_name);
		for (size_t i = 0; i < results.size(); i++) {
			os << results[i] << endl;
		}
		os.close();
	}

	void WriteResults(const string &folder_name) {
		mkdir(folder_name.c_str(),
				S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
		WriteResultsToFile(perfect_matches_, folder_name + "/perfect_matches.txt");
		WriteResultsToFile(good_matches_, folder_name + "/good_matches.txt");
		WriteResultsToFile(mismatches_, folder_name + "/mismatches.txt");
		WriteResultsToFile(imperfect_matches_,
				folder_name + "/imperfect_matches.txt");
	}
};


}

#endif /* OMNI_TOOLS_HPP_ */

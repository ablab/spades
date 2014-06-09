//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNI_TOOLS_HPP_
#define OMNI_TOOLS_HPP_

#include "omni_utils.hpp"
#include "simple_tools.hpp"

#include "path_helper.hpp"

#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

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
	bool safe_merging_;

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

//do not use without checks:)
	EdgeId CompressWithoutChecks(VertexId v) {

		EdgeId e = graph_.GetUniqueOutgoingEdge(v);
		EdgeId start_edge = e;
		while (GoUniqueWayBackward(e) && e != start_edge
				&& !graph_.RelatedVertices(graph_.EdgeStart(e),
						graph_.EdgeEnd(e))) {
		}
		vector<EdgeId> mergeList;
		//		e = graph_.conjugate(e);
		start_edge = e;
		do {
			mergeList.push_back(e);
		} while (GoUniqueWayForward(e) && e != start_edge
				&& !graph_.RelatedVertices(graph_.EdgeStart(e),
						graph_.EdgeEnd(e)));
		EdgeId new_edge = graph_.MergePath(mergeList, safe_merging_);
		TRACE("Vertex compressed and is now part of edge "
						<< graph_.str(new_edge));
		return new_edge;

	}

	//todo use graph method!
	bool CanCompressVertex(VertexId v) {
		if (!graph_.CheckUniqueOutgoingEdge(v)
			|| !graph_.CheckUniqueIncomingEdge(v)) {
			TRACE(
					"Vertex "
							<< graph_.str(v)
							<< " judged NOT compressible. Proceeding to the next vertex");
			TRACE("Processing vertex " << graph_.str(v) << " finished");
			return false;
		}
		return true;
	}
public:
	Compressor(Graph &graph, bool safe_merging = true) :
			graph_(graph),
			safe_merging_(safe_merging) {
	}

	/**
	 * Method compresses longest possible path, containing given vertex.
	 * @param vertex to be compressed as part of a path
	 * @return true if vertex can be compressed and false otherwise
	 */
	bool CompressVertex(VertexId v) {
		TRACE("Processing vertex " << graph_.str(v) << " started");
		if (! CanCompressVertex(v)) {
			return false;
		}
		TRACE("Vertex " << graph_.str(v) << " judged compressible");
		CompressWithoutChecks(v);
		return true;
	}

	EdgeId CompressVertexEdgeId(VertexId v){
		TRACE("Processing vertex " << graph_.str(v) << " started");
		if (! CanCompressVertex(v)) {
			return EdgeId(0);
		}
		TRACE("Vertex " << graph_.str(v) << " judged compressible");
		return CompressWithoutChecks(v);

	}

	/**
	 * Method compresses all vertices which can be compressed.
	 */
	void CompressAllVertices() {
		TRACE("Vertex compressing started");
		//SmartVertexIterator<Graph> end = graph_.SmartVertexEnd();

		//todo is it still true???
		//in current implementation will work incorrectly if smart iterator won't give vertex and its conjugate
		//(in case of self-conjugate edges)
		vector<VertexId> to_compress;
		for (VertexId v : graph_) {
			if (CanCompressVertex(v)) {
				to_compress.push_back(v);
			}
		}
		for (auto it = SmartSetIterator<Graph, VertexId>(graph_, to_compress.begin(), to_compress.end()); !it.IsEnd(); ++it) {
			CompressVertex(*it);
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
	DECL_LOGGER("Cleaner")
};

template<class Graph>
class IsolatedEdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	size_t max_length_;
	double max_coverage_;
	size_t max_length_any_cov_;

	bool IsTerminalVertex(VertexId v) {
		return g_.IncomingEdgeCount(v) + g_.OutgoingEdgeCount(v) == 1;
	}

public:
	IsolatedEdgeRemover(Graph& g, size_t max_length, double max_coverage, size_t max_length_any_cov) :
			g_(g), max_length_(max_length), max_coverage_(max_coverage), max_length_any_cov_(max_length_any_cov) {
	}

	size_t RemoveIsolatedEdges() {
        size_t cnt = 0;
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (IsTerminalVertex(g_.EdgeStart(*it))
					&& IsTerminalVertex(g_.EdgeEnd(*it))
					&& ((g_.length(*it) <= max_length_ && g_.coverage(*it) <= max_coverage_)
					|| g_.length(*it) <=max_length_any_cov_)) {
				g_.DeleteEdge(*it);
                cnt += 1;

			}
		}
		Cleaner<Graph> cleaner(g_);
		cleaner.Clean();

        return cnt;
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
			DEBUG("Added vertex " << new_vertex);
			VerticesCopies.insert(make_pair(*iter, new_vertex));
		}
		for (auto iter = graph_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
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
		StructureCopy<CopyGraph>(new_graph, VerticesCopies, EdgesCopies);
	}
	template<class CopyGraph>
	void TotalCopy(CopyGraph &new_graph) {
		typedef typename CopyGraph::EdgeId NewEdgeId;
		typedef typename CopyGraph::VertexId NewVertexId;
		map<EdgeId, NewEdgeId> EdgesCopies;
		map<VertexId, NewVertexId> VerticesCopies;
		StructureCopy<CopyGraph>(new_graph, VerticesCopies, EdgesCopies);
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
		for (auto iter = graph_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
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
	TrivialEdgePairChecker(const Graph &graph, size_t bound = (size_t) - 1) :
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
class AvgCovereageCounter {
private:
	const Graph &graph_;
	const size_t min_length_;
public:
	AvgCovereageCounter(const Graph &graph, size_t min_length = 0) :
			graph_(graph), min_length_(min_length) {
	}

	double Count() const {
		double cov = 0;
		size_t length = 0;
		for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
			if (graph_.length(*it) >= min_length_) {
				cov += graph_.coverage(*it) * (double) graph_.length(*it);
				length += graph_.length(*it);
			}
		}
		if (length == 0)
			return 0.;
		return cov / (double) length;
	}
};

template<class Graph>
class ErroneousConnectionThresholdFinder {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    size_t backet_width_;

    bool IsInteresting(EdgeId e) const {
        if (graph_.length(e) > graph_.k() + 1)
            return false;

        if (graph_.OutgoingEdgeCount(graph_.EdgeStart(e)) < 2 ||
            graph_.IncomingEdgeCount(graph_.EdgeEnd(e)) < 2)
            return false;

        std::vector<EdgeId> v1;
        push_back_all(v1, graph_.OutgoingEdges(graph_.EdgeStart(e)));
        std::vector<EdgeId> v2;
        push_back_all(v2, graph_.IncomingEdges(graph_.EdgeEnd(e)));
        bool eq = (v1.size() == 2 && v2.size() == 2) && ((v1[0] == v2[0] && v1[1] == v2[1])	|| (v1[0] == v2[1] && v1[0] == v2[1]));
        return !eq;
    }

    double weight(size_t value, const map<size_t, size_t> &histogram,
                  size_t backet_width) const {
        double result = 0;
        for (size_t i = 0; i < backet_width && value + i < histogram.size(); i++) {
            result += (double) (getValue(value + i, histogram) * std::min(i + 1, backet_width - i));
        }
        return result;
    }

    double Median(double thr = 500.0) const {
        vector<double> coverages;
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (graph_.length(*it) > thr)
                coverages.push_back(graph_.coverage(*it));
        }

        auto middle_it = coverages.begin() + coverages.size() / 2;
#ifdef USE_GLIBCXX_PARALLEL
        __gnu_parallel::nth_element(coverages.begin(), middle_it, coverages.end());
#else
        std::nth_element(coverages.begin(), middle_it, coverages.end());
#endif
        return coverages[coverages.size() / 2];
    }

    size_t getValue(size_t arg, const map<size_t, size_t> &ssmap) const {
        auto it = ssmap.find(arg);
        if (it == ssmap.end())
            return 0;
        else
            return it->second;
    }

public:
    ErroneousConnectionThresholdFinder(const Graph &graph, size_t backet_width = 0) :
            graph_(graph), backet_width_(backet_width) {
    }

    double AvgCoverage() const {
        double cov = 0;
        double length = 0;
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            cov += graph_.coverage(*it) * (double) graph_.length(*it);
            length += (double) graph_.length(*it);
        }
        return cov / length;
    }

    std::map<size_t, size_t> ConstructHistogram() const {
        std::map<size_t, size_t> result;
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (IsInteresting(*it))
                result[(size_t)graph_.coverage(*it)]++;
        }
        return result;
    }

    double FindThreshold(const map<size_t, size_t> &histogram) const {
        size_t backet_width = backet_width_;
        if (backet_width == 0) {
            backet_width = (size_t)(0.3 * AvgCovereageCounter<Graph>(graph_).Count() + 5);
        }
        size_t size = 0;
        if (histogram.size() != 0)
            size = histogram.rbegin()->first + 1;
        INFO("Bucket size: " << backet_width);
        size_t cnt = 0;
        for (size_t i = 1; i + backet_width < size; i++) {
            if (weight(i, histogram, backet_width) > weight(i - 1, histogram, backet_width))
                cnt++;

            if (i > backet_width &&
                weight(i - backet_width,     histogram, backet_width) >
                weight(i - backet_width - 1, histogram, backet_width)) {
                cnt--;
            }
            if (2 * cnt >= backet_width)
                return (double) i;

        }
        INFO("Proper threshold was not found. Threshold set to 0.1 of average coverage");
        return 0.1 * AvgCovereageCounter<Graph>(graph_).Count();
    }

    double FindThreshold() const {
        INFO("Finding threshold started");
        std::map<size_t, size_t> histogram = ConstructHistogram(/*weights*/);
        for (size_t i = 0; i < histogram.size(); i++) {
            TRACE(i << " " << histogram[i]);
        }
        double result = FindThreshold(histogram);
        INFO("Average edge coverage: " << AvgCoverage());
        INFO("Graph threshold: " << result);
        result = std::max(AvgCoverage(), result);
        INFO("Threshold finding finished. Threshold is set to " << result);
        return result;
    }
private:
    DECL_LOGGER("ThresholdFinder");
};

}

#endif /* OMNI_TOOLS_HPP_ */

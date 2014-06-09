//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard_base.hpp"
#include "dijkstra_tools/dijkstra_helper.hpp"

namespace omnigraph {

template<class Graph>
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	std::stringstream ss;
	for (size_t i = 0; i < edges.size(); ++i) {
		ss << delim << g.str(edges[i]) << " (" << g.length(edges[i]) << ")";
		delim = " -> ";
	}
	return ss.str();
}

template<class Graph>
class PathProcessor {

	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<EdgeId> Path;
	typedef typename DijkstraHelper<Graph>::BoundedDijkstra DijkstraT;

public:
	class Callback {

	public:
		virtual ~Callback() {
		}

		virtual void Flush() {};

		virtual void HandleReversedPath(
				const vector<EdgeId>& reversed_path) = 0;

	protected:
		Path ReversePath(const Path& path) const {
			Path result;
			for (auto I = path.rbegin(), E = path.rend(); I != E; ++I)
				result.push_back(*I);
			return result;
		}
	};

	int Go(VertexId v, const size_t min_len, size_t cur_len,
			const DijkstraT& dsts_to_start) {
		int error_code = 0;
		//TRACE("Now in the vertex " << g_.int_id(v));
		if (++call_cnt_ >= MAX_CALL_CNT) {
			//DEBUG("Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
			return 1;
		}

		if (!dsts_to_start.DistanceCounted(v)
				|| dsts_to_start.GetDistance(v) + cur_len > max_len_) {
			//if (!dsts_to_start.DistanceCounted(v)) {
			//TRACE("Distance not counted yet");
			//}
			//else
			//TRACE("Shortest distance from this vertex is " << dsts_to_start.GetDistance(v)
			//<< " and sum with current path length " << cur_len
			//<< " exceeded max length " << max_len_);
			return 0;
		}

		if (v == start_ && cur_len >= min_len) {
			//TRACE("New path found: " << PrintPath(g_, path_));
			callback_->HandleReversedPath(path_);
		}
		TRACE("Iterating through incoming edges of vertex " << g_.int_id(v))
		//TODO: doesn`t work with parallel simplification
		for (auto I = g_.in_begin(v), E = g_.in_end(v); I != E; ++I) {
			EdgeId edge = *I;
//		BOOST_FOREACH(EdgeId edge, g_.IncomingEdges(v)) {
			path_.push_back(edge);
			error_code |= Go(g_.EdgeStart(edge), min_len,
					cur_len + g_.length(edge), dsts_to_start);
			path_.pop_back();
		}
		//TRACE("Processing vertex " << g_.int_id(v) << " finished");
		return error_code;
	}

	// constructor for paths between start vertex and a set of @end_points
	PathProcessor(const Graph& g, vector<size_t> min_lens, size_t max_len,
			VertexId start, vector<VertexId> end_points, Callback& callback) :
			g_(g), min_lens_(min_lens), max_len_(max_len), start_(start), end_points_(
					end_points), dijkstra_(
					DijkstraHelper<Graph>::CreateBoundedDijkstra(g, max_len, MAX_DIJKSTRA_VERTICES)),
					callback_(&callback), call_cnt_(0) {
		dijkstra_.run(start);
	}

	// constructor when we have only one @end_point
	PathProcessor(const Graph& g, size_t min_len, size_t max_len,
			VertexId start, VertexId end_point, Callback& callback) :
			g_(g), max_len_(max_len), start_(start), dijkstra_(
					DijkstraHelper<Graph>::CreateBoundedDijkstra(g, max_len, MAX_DIJKSTRA_VERTICES)),
					callback_(&callback), call_cnt_(0) {
		min_lens_.push_back(min_len);
		end_points_.push_back(end_point);
		dijkstra_.run(start);
	}

	// dfs from the end vertices
	int Process() {
		int error_code = 0;

		if (dijkstra_.VertexLimitExceeded()) {
			TRACE("dijkstra : vertex limit exceeded");
			error_code = 2;
		}

		TRACE("Start vertex is " << g_.int_id(start_));
		for (size_t i = 0; i < end_points_.size(); ++i) {
			call_cnt_ = 0;
			VertexId current_end_ = end_points_[i];
			TRACE("Bounds are " << min_lens_[i] << " " << max_len_);
			TRACE("Current end vertex " << g_.int_id(current_end_));
			error_code |= Go(current_end_, min_lens_[i], 0, dijkstra_);
			callback_->Flush();
		}
		return error_code; // 3 two mistakes, 2 bad dijkstra, 1 bad dfs, 0 = okay
	}

	//todo remove setters
	void SetMinLens(const vector<size_t>& new_min_lens) {
		min_lens_ = new_min_lens;
	}

	void SetMinLens(vector<size_t> && new_min_lens) {
		min_lens_ = new_min_lens;
	}

	void SetMaxLen(size_t new_max_len) {
		max_len_ = new_max_len;
	}

	void SetEndPoints(const vector<VertexId>& new_end_points) {
		end_points_ = new_end_points;
	}

	void SetEndPoints(vector<VertexId>&& new_end_points) {
		end_points_ = new_end_points;
	}

	void SetCallback(Callback* new_callback) {
		callback_ = new_callback;
	}

	void ResetCallCount() {
		call_cnt_ = 0;
	}

private:
	static const size_t MAX_CALL_CNT = 3000;
	static const size_t MAX_DIJKSTRA_VERTICES = 3000;

	const Graph& g_;
	vector<size_t> min_lens_;
	size_t max_len_;
	VertexId start_;
	vector<VertexId> end_points_;
	DijkstraT dijkstra_;
	Callback* callback_;
	Path path_;
	size_t call_cnt_;

	DECL_LOGGER("PathProcessor")
};

template<class Graph>
class CompositeCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<EdgeId> Path;

public:
	void AddProcessor(typename PathProcessor<Graph>::Callback& processor) {
		processors_.push_back(&processor);
	}

	virtual void Flush() {
		for (auto it = processors_.begin(); it != processors_.end(); ++it) {
			(*it)->Flush();
		}
	}

	virtual void HandleReversedPath(const Path& path) {
		for (auto it = processors_.begin(); it != processors_.end(); ++it) {
			(*it)->HandleReversedPath(path);
		}
	}

private:
	vector<typename PathProcessor<Graph>::Callback*> processors_;
};

template<class Graph>
class PathStorageCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<EdgeId> Path;

public:
	PathStorageCallback(const Graph& g) :
			g_(g) {
	}

	virtual void Flush() {
		all_paths_.push_back(cur_paths_);
		cur_paths_.clear();
	}

	virtual void HandleReversedPath(const vector<EdgeId>& path) {
		cur_paths_.push_back(this->ReversePath(path));
	}

	size_t size(size_t k = 0) const {
		return all_paths_[k].size();
	}

	const vector<Path>& paths(size_t k = 0) const {
		return all_paths_[k];
	}

private:
	const Graph& g_;
	vector<vector<Path>> all_paths_;
	vector<Path> cur_paths_;
};

template<class Graph>
class NonEmptyPathCounter: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<EdgeId> Path;

public:
	NonEmptyPathCounter(const Graph& g) :
			g_(g), count_(0) {
	}

	virtual void Flush() {
		all_paths_.push_back(cur_paths_);
		counts_.push_back(count_);
		cur_paths_.clear();
	}

	virtual void HandleReversedPath(const Path& path) {
		if (path.size() > 0) {
			++count_;
			cur_paths_.push_back(this->ReversePath(path));
		}
	}

	size_t count(size_t k = 0) const {
		return counts_[k];
	}

	const vector<Path>& paths(size_t k = 0) const {
		return all_paths_[k];
	}

private:
	const Graph& g_;
	vector<size_t> counts_;
	size_t count_;
	vector<vector<Path> > all_paths_;
	vector<Path> cur_paths_;
};

template<class Graph>
class VertexLabelerCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<EdgeId> Path;

public:
	VertexLabelerCallback(const Graph& g) :
			g_(g), count_(0) {
	}

	virtual void Flush() {
		all_vertices_.push_back(vertices_);
		vertices_.clear();
		counts_.push_back(count_);
	}

	virtual void HandleReversedPath(const Path& path) {
		for (auto it = path.rbegin(); it != path.rend(); ++it) {
			if (path.size() > 0) {
				vertices_.insert(g_.EdgeStart(*it));
				vertices_.insert(g_.EdgeEnd(*it));
				++count_;
			}
		}
	}

	const set<VertexId>& vertices(size_t k = 0) const {
		return all_vertices_[k];
	}

	size_t count(size_t k = 0) const {
		return counts_[k];
	}

private:
	Graph& g_;
	vector<size_t> counts_;
	vector<set<VertexId>> all_vertices_;
	size_t count_;
	set<VertexId> vertices_;
};

template<class Graph>
class DistancesLengthsCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef vector<EdgeId> Path;

public:
	DistancesLengthsCallback(const Graph& g) :
			g_(g) {
	}

	virtual void Flush() {
		all_distances_.push_back(distances_);
		distances_.clear();
	}

	virtual void HandleReversedPath(const Path& path) {
		size_t path_length = PathLength(path);
		distances_.insert(path_length);
	}

	vector<size_t> distances(size_t k = 0) const {
		const set<size_t>& tmp = all_distances_[k];
		return vector<size_t>(tmp.begin(), tmp.end());
	}

private:
	size_t PathLength(const Path& path) const {
		size_t res = 0;
		for (auto I = path.begin(); I != path.end(); ++I)
			res += g_.length(*I);
		return res;
	}

	const Graph& g_;
	set<size_t> distances_;
	vector<set<size_t>> all_distances_;

	DECL_LOGGER("DistancesLengthsCallback");
};

template<class Graph>
class MappingPathFixer {
public:

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    MappingPathFixer(const Graph& graph)
            : g_(graph) {
    }

    bool CheckContiguous(const vector<typename Graph::EdgeId>& path) const {
        for (size_t i = 1; i < path.size(); ++i) {
            if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i]))
                return false;
        }
        return true;
    }

    //todo seems that we don't need optional here any more
    vector<EdgeId> TryFixPath(const vector<EdgeId>& edges) const {
        vector<EdgeId> answer;
        if (edges.empty()) {
            //          WARN("Mapping path was empty");
            return vector<EdgeId>();
        }
        //      VERIFY(edges.size() > 0);
        answer.push_back(edges[0]);
        for (size_t i = 1; i < edges.size(); ++i) {
            if (g_.EdgeEnd(edges[i - 1]) != g_.EdgeStart(edges[i])) {
                vector<EdgeId> closure = TryCloseGap(g_.EdgeEnd(edges[i - 1]),
                                                     g_.EdgeStart(edges[i]));
                answer.insert(answer.end(), closure.begin(), closure.end());
                //                  make_dir("assembly_compare/tmp");
                //                  WriteComponentsAroundEdge(graph_,
                //                          graph_.IncomingEdges(v1).front(),
                //                          "assembly_compare/tmp/failed_close_gap_from.dot",
                //                          *DefaultColorer(graph_),
                //                          LengthIdGraphLabeler<Graph>(g_));
                //                  WriteComponentsAroundEdge(graph_,
                //                          g_.OutgoingEdges(v2).front(),
                //                          "assembly_compare/tmp/failed_close_gap_to.dot",
                //                          *DefaultColorer(g_),
                //                          LengthIdGraphLabeler<Graph>(g_));
                //                  VERIFY(false);
            }
            answer.push_back(edges[i]);
        }
        return answer;
    }

    vector<EdgeId> DeleteSameEdges(const vector<EdgeId>& path) const {
        vector<EdgeId> result;
        if (path.empty()) {
            return result;
        }
        result.push_back(path[0]);
        for (size_t i = 1; i < path.size(); ++i) {
            if (path[i] != result[result.size() - 1]) {
                result.push_back(path[i]);
            }
        }
        return result;
    }

private:
    vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2) const {
        if (v1 == v2)
            return vector<EdgeId>();
        TRACE(
                "Trying to close gap between v1=" << g_.int_id(v1) << " and v2=" << g_.int_id(v2));
        PathStorageCallback<Graph> path_store(g_);
        //todo reduce value after investigation
        PathProcessor<Graph> path_processor(g_, 0, 70, v1, v2, path_store);
        path_processor.Process();

        if (path_store.size() == 0) {
            TRACE("Failed to find closing path");
            //          TRACE("Failed to close gap between v1=" << graph_.int_id(v1)
            //                          << " (conjugate "
            //                          << graph_.int_id(g_.conjugate(v1))
            //                          << ") and v2=" << g_.int_id(v2)
            //                          << " (conjugate "
            //                          << g_.int_id(g_.conjugate(v2)) << ")");
            //          return boost::none;
            return vector<EdgeId>();
        } else if (path_store.size() == 1) {
            TRACE("Unique closing path found");
        } else {
            TRACE("Several closing paths found, first chosen");
        }
        vector<EdgeId> answer = path_store.paths().front();
        TRACE("Gap closed");
        TRACE( "Cumulative closure length is " << CumulativeLength(g_, answer));
        return answer;
    }
    const Graph& g_;
};



}

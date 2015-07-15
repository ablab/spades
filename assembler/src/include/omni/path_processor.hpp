//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "standard_base.hpp"
#include "adt/bag.hpp"
#include "dijkstra_tools/dijkstra_helper.hpp"

namespace omnigraph {

template<class Graph>
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	std::stringstream ss;
	for (size_t i = 0; i < edges.size(); ++i) {
		ss << delim << g.str(edges[i]);
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

    void Push(EdgeId e, VertexId start_v) {
        TRACE("Pushing edge " << g_.str(e));
        curr_len_ += g_.length(e);
        reversed_edge_path_.push_back(e);
        vertex_cnts_.put(start_v);
    }

    void Pop() {
        VERIFY(!reversed_edge_path_.empty());
        EdgeId e = reversed_edge_path_.back();
        size_t len = g_.length(e);
        VERIFY(curr_len_ >= len);

        TRACE("Popping edge " << g_.str(e));
        vertex_cnts_.take(g_.EdgeStart(e));
        reversed_edge_path_.pop_back();
        curr_len_ -= len;
    }

    bool CanGo(EdgeId e, VertexId start_v) {
//        VertexId v = g_.EdgeStart(e);
        //if (!dijkstra_.DistanceCounted(v)) {
        //TRACE("Distance not counted yet");
        //}
        //else
        //TRACE("Shortest distance from this vertex is " << dijkstra_.GetDistance(v)
        //<< " and sum with current path length " << cur_len
        //<< " exceeded max length " << max_len_);
        if (!dijkstra_.DistanceCounted(start_v))
            return false;
        if (dijkstra_.GetDistance(start_v) + g_.length(e) + curr_len_ > max_len_)
            return false;
        if (vertex_cnts_.mult(start_v) >= MAX_VERTEX_USAGE)
            return false;
        return true;
    }

    //returns true iff limits were exceeded
    bool Go(VertexId v, const size_t min_len) {
        TRACE("Got to vertex " << g_.str(v));
        if (++call_cnt_ >= MAX_CALL_CNT) {
            TRACE("Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
            return true;
        }

        if (v == start_ && curr_len_ >= min_len) {
            //TRACE("New path found: " << PrintPath(g_, path_));
            callback_->HandleReversedPath(reversed_edge_path_);
        }

        TRACE("Iterating through incoming edges of vertex " << g_.int_id(v))
        //TODO: doesn`t work with parallel simplification
        vector<EdgeId> incoming;
        incoming.reserve(4);
        std::copy_if(g_.in_begin(v), g_.in_end(v), std::back_inserter(incoming), [&] (EdgeId e) {
            return dijkstra_.DistanceCounted(g_.EdgeStart(e));
        });

        std::sort(incoming.begin(), incoming.end(), [&] (EdgeId e1, EdgeId e2) {
            return dijkstra_.GetDistance(g_.EdgeStart(e1)) < dijkstra_.GetDistance(g_.EdgeStart(e2));
        });

        for (EdgeId e : incoming) {
            VertexId start_v = g_.EdgeStart(e);
            if (CanGo(e, start_v)) {
                Push(e, start_v);
                bool exceeded_limits = Go(start_v, min_len);
                Pop();
                if (exceeded_limits)
                    return true;
            }
        }
        //TRACE("Processing vertex " << g_.int_id(v) << " finished");
        return false;
    }

public:
    class Callback {

    public:
        virtual ~Callback() {
        }

        virtual void Flush() {
        }

        virtual void HandleReversedPath(const vector<EdgeId>& reversed_path) = 0;


    protected:
        Path ReversePath(const Path& path) const {
            Path result;
            for (auto I = path.rbegin(), E = path.rend(); I != E; ++I)
                result.push_back(*I);
            return result;
        }
    };

    // constructor for paths between start vertex and a set of @end_points
    PathProcessor(const Graph& g, const vector<size_t>& min_lens, size_t max_len, VertexId start,
                  const vector<VertexId>& end_points, Callback& callback)
            : g_(g),
              min_lens_(min_lens),
              max_len_(max_len),
              start_(start),
              end_points_(end_points),
              dijkstra_(DijkstraHelper<Graph>::CreateBoundedDijkstra(g, max_len, MAX_DIJKSTRA_VERTICES)),
              callback_(&callback),
              curr_len_(0),
              call_cnt_(0) {
        TRACE("Dijkstra launched");
        dijkstra_.run(start);
        reversed_edge_path_.reserve(MAX_CALL_CNT);
        TRACE("Dijkstra finished");
    }

    // constructor when we have only one @end_point
    PathProcessor(const Graph& g, size_t min_len, size_t max_len, VertexId start, VertexId end_point, Callback& callback)
            : g_(g),
              max_len_(max_len),
              start_(start),
              dijkstra_(DijkstraHelper<Graph>::CreateBoundedDijkstra(g, max_len, MAX_DIJKSTRA_VERTICES)),
              callback_(&callback),
              curr_len_(0),
              call_cnt_(0) {
        TRACE("Dijkstra launched");
        min_lens_.push_back(min_len);
        end_points_.push_back(end_point);
        dijkstra_.run(start);
        reversed_edge_path_.reserve(MAX_CALL_CNT);
        TRACE("Dijkstra finished");
    }

    // dfs from the end vertices
    // 3 two mistakes, 2 bad dijkstra, 1 some bad dfs, 0 = okay
    int Process() {
        TRACE("Process launched");
        int error_code = 0;

        if (dijkstra_.VertexLimitExceeded()) {
            TRACE("dijkstra : vertex limit exceeded");
            error_code = 2;
        }

        TRACE("Start vertex is " << g_.int_id(start_));
        for (size_t i = 0; i < end_points_.size(); ++i) {
            VERIFY(curr_len_ == 0);
            VERIFY(vertex_cnts_.size() == 0);
            call_cnt_ = 0;
            VertexId current_end = end_points_[i];
            TRACE("Bounds are " << min_lens_[i] << " " << max_len_);
            TRACE("Current end vertex " << g_.int_id(current_end));
            vertex_cnts_.put(current_end);
            error_code |= int(Go(current_end, min_lens_[i]));
            vertex_cnts_.take(current_end);
            callback_->Flush();
        }
        TRACE("Process finished with error code " << error_code);
        return error_code;
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

    void SetEndPoints(vector<VertexId> && new_end_points) {
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
    static const size_t MAX_VERTEX_USAGE = 5;

    const Graph& g_;
    vector<size_t> min_lens_;
    size_t max_len_;
    VertexId start_;
    vector<VertexId> end_points_;
    DijkstraT dijkstra_;
    Callback* callback_;

    Path reversed_edge_path_;
    bag<VertexId> vertex_cnts_;
    size_t curr_len_;
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

template<class Graph, class Comparator>
class BestPathStorage: public PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef vector<EdgeId> Path;
public:
    BestPathStorage(const Graph& g, Comparator comparator) :
            g_(g), cnt_(0), comparator_(comparator) {
    }

    virtual void HandleReversedPath(const vector<EdgeId>& path) {
        cnt_++;
        if(best_path_.size() == 0 || comparator_(path, best_path_))
            best_path_ = path;
    }

    vector<EdgeId> BestPath() const {
        return best_path_;
    }

    size_t size() const {
        return cnt_;
    }

private:
    const Graph& g_;
    size_t cnt_;
    Comparator comparator_;
    vector<vector<Path>> best_path_;
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

    Path<EdgeId> TryFixPath(const Path<EdgeId>& path, size_t length_bound = 70) const {
        return Path<EdgeId>(TryFixPath(path.sequence(), length_bound), path.start_pos(), path.end_pos());
    }

    vector<EdgeId> TryFixPath(const vector<EdgeId>& edges, size_t length_bound = 70) const {
        vector<EdgeId> answer;
        if (edges.empty()) {
            //          WARN("Mapping path was empty");
            return vector<EdgeId>();
        }
        answer.push_back(edges[0]);
        for (size_t i = 1; i < edges.size(); ++i) {
            if (g_.EdgeEnd(edges[i - 1]) != g_.EdgeStart(edges[i])) {
                vector<EdgeId> closure = TryCloseGap(g_.EdgeEnd(edges[i - 1]),
                                                     g_.EdgeStart(edges[i]),
                                                     length_bound);
                answer.insert(answer.end(), closure.begin(), closure.end());
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
    vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2, size_t length_bound) const {
        if (v1 == v2)
            return vector<EdgeId>();
        TRACE(
                "Trying to close gap between v1=" << g_.int_id(v1) << " and v2=" << g_.int_id(v2));
        PathStorageCallback<Graph> path_store(g_);
        //todo reduce value after investigation
        PathProcessor<Graph> path_processor(g_, 0, length_bound, v1, v2, path_store);
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

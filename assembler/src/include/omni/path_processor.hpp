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
const string PrintPath(const Graph& g, const vector<typename Graph::EdgeId>& edges) {
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

private:

    class Traversal {
        const PathProcessor& outer_;
        VertexId end_;
        size_t min_len_;
        size_t max_len_;
        Callback& callback_;
        size_t edge_depth_bound_;

        size_t curr_len_;
        size_t curr_depth_;
        size_t call_cnt_;
        Path reversed_edge_path_;
        bag<VertexId> vertex_cnts_;

        const Graph& g_;
        const DijkstraT& dijkstra_;

        void Push(EdgeId e, VertexId start_v) {
            TRACE("Pushing edge " << g_.str(e));
            curr_len_ += g_.length(e);
            curr_depth_++;
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
            curr_depth_--;
        }

        bool CanGo(EdgeId e, VertexId start_v) {
            if (!dijkstra_.DistanceCounted(start_v))
                return false;
            if (dijkstra_.GetDistance(start_v) + g_.length(e) + curr_len_ > max_len_)
                return false;
            if (curr_depth_ >= edge_depth_bound_)
                return false;
            if (vertex_cnts_.mult(start_v) >= PathProcessor::MAX_VERTEX_USAGE)
                return false;
            return true;
        }

        bool Go(VertexId v, const size_t min_len) {
            TRACE("Got to vertex " << g_.str(v));
            if (++call_cnt_ >= PathProcessor::MAX_CALL_CNT) {
                TRACE("Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
                return true;
            }

            if (v == outer_.start_ && curr_len_ >= min_len) {
                //TRACE("New path found: " << PrintPath(g_, path_));
                callback_.HandleReversedPath(reversed_edge_path_);
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
            return false;
        }

    public:
        Traversal(const PathProcessor& outer, VertexId end,
                  size_t min_len, size_t max_len,
                  Callback& callback, size_t edge_depth_bound) :
            outer_(outer), end_(end),
            min_len_(min_len), max_len_(max_len),
            callback_(callback),
            edge_depth_bound_(edge_depth_bound),
            curr_len_(0), curr_depth_(0), call_cnt_(0), 
            g_(outer.g_),
            dijkstra_(outer.dijkstra_) {
            reversed_edge_path_.reserve(PathProcessor::MAX_CALL_CNT);
            vertex_cnts_.put(end_);
        }

        //returns true iff limits were exceeded
        bool Go() {
            if (!dijkstra_.DistanceCounted(end_) || dijkstra_.GetDistance(end_) > max_len_) {
                return false;
            }

            bool code = Go(end_, min_len_);
            VERIFY(curr_len_ == 0);
            VERIFY(curr_depth_ == 0);
            vertex_cnts_.take(end_);
            VERIFY(vertex_cnts_.size() == 0);
            return code;
        }
    };

    friend class Traversal;

public:

    PathProcessor(const Graph& g, VertexId start, size_t length_bound) :
    		  g_(g),
              start_(start),
              dijkstra_(DijkstraHelper<Graph>::CreateBoundedDijkstra(g, length_bound, MAX_DIJKSTRA_VERTICES)) {
        TRACE("Dijkstra launched");
        dijkstra_.Run(start);
        TRACE("Dijkstra finished");
    }

    // dfs from the end vertices
    // 3 two mistakes, 2 bad dijkstra, 1 some bad dfs, 0 = okay
    int Process(VertexId end, size_t min_len, size_t max_len, Callback& callback, size_t edge_depth_bound = -1ul) const {
        TRACE("Process launched");
        int error_code = 0;

        if (dijkstra_.VertexLimitExceeded()) {
            TRACE("dijkstra : vertex limit exceeded");
            error_code = 2;
        }

        TRACE("Start vertex is " << g_.str(start_));
        TRACE("Bounds are " << min_len << " " << max_len);
        TRACE("End vertex " << g_.str(end));

        Traversal traversal(*this, end, min_len, max_len, callback, edge_depth_bound);
        error_code |= int(traversal.Go());

        callback.Flush();
        TRACE("Process finished with error code " << error_code);
        return error_code;
    }

private:
    static const size_t MAX_CALL_CNT = 3000;
    static const size_t MAX_DIJKSTRA_VERTICES = 3000;
    static const size_t MAX_VERTEX_USAGE = 5;

    const Graph& g_;
    VertexId start_;
    DijkstraT dijkstra_;

    DECL_LOGGER("PathProcessor")
};

template<class Graph>
int ProcessPaths(const Graph& g, size_t min_len, size_t max_len,
                 typename Graph::VertexId start, typename Graph::VertexId end,
                 typename PathProcessor<Graph>::Callback& callback, size_t max_edge_cnt = -1ul) {
    PathProcessor<Graph> processor(g, start, max_len);
    return processor.Process(end, min_len, max_len, callback, max_edge_cnt);
}

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
	    VERIFY(k < all_distances_.size());
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

}

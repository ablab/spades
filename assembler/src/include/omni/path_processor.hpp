#pragma once

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

 public:
  class Callback {

   public:
    virtual ~Callback()
    {
    }

    virtual void Flush() = 0;

    virtual void HandleReversedPath(const vector<EdgeId>& reversed_path) = 0;
   
   protected:
    Path ReversePath(const Path& path) const {
      Path result;
      for (auto I = path.rbegin(), E = path.rend(); I != E; ++I)
        result.push_back(*I);
      return result;
    }
  };

  int Go(VertexId v, const size_t min_len, size_t cur_len, Dijkstra<Graph>& dsts_to_start)
  {
    int error_code = 0;
    TRACE("Now in the vertex " << g_.int_id(v));
    if (++call_cnt_ >= MAX_CALL_CNT) {
      DEBUG("Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
      return 1;
    }

    if (!dsts_to_start.DistanceCounted(v) ||
         dsts_to_start.GetDistance(v) + cur_len > max_len_)
    {
      if (!dsts_to_start.DistanceCounted(v)) {
        TRACE("Distance not counted yet");
      }
      else
        TRACE("Shortest distance from this vertex is " << dsts_to_start.GetDistance(v)
            << " and sum with current path length " << cur_len
            << " exceeded max length " << max_len_);
      return 0;
    }

    if (v == start_ && cur_len >= min_len) {
      TRACE("New path found: " << PrintPath(g_, path_));
      callback_.HandleReversedPath(path_);
    }
    TRACE("Iterating through incoming edges of vertex " << g_.int_id(v))
    //TODO: doesn`t work with parallel simprification
    //for (auto I = g_.in_begin(v), E = g_.in_end(v); I != E; ++I) {
      //EdgeId edge = *I;
    BOOST_FOREACH(EdgeId edge, g_.IncomingEdges(v)) {
      path_.push_back(edge);
      error_code |= Go(g_.EdgeStart(edge), min_len, cur_len + g_.length(edge),  dsts_to_start);
      path_.pop_back();
    }
    TRACE("Processing vertex " << g_.int_id(v) << " finished");
    return error_code;
  }

  // constructor for paths between start vertex and a set of @end_points
  PathProcessor(const Graph& g, const vector<size_t>& min_lens, size_t max_len,
                 VertexId start, const vector<VertexId>& end_points, Callback& callback) :
                 g_(g), min_lens_(min_lens), max_len_(max_len),
                 start_(start), end_points_(end_points), callback_(callback), call_cnt_(0)
  {
  }

  // constructor when we have only one @end_point
  PathProcessor(const Graph& g, size_t min_len, size_t max_len,
                 VertexId start, VertexId end_point, Callback& callback) :
                 g_(g), min_lens_({min_len}), max_len_(max_len),
                 start_(start), end_points_({end_point}), callback_(callback), call_cnt_(0)
  {
  }

  // dfs from the end vertices
   int Process() {
    int error_code = 0;
    ReliableBoundedDijkstra<Graph> dijkstra(g_, max_len_, MAX_DIJKSTRA_VERTICES);
    dijkstra.run(start_);

    if (dijkstra.VertexLimitExceeded()) {
      TRACE("dijkstra : vertex limit exceeded");
      error_code = 2;
    }

    TRACE("Start vertex is " << g_.int_id(start_));
    for (size_t i = 0; i < end_points_.size(); ++i) {
      VertexId current_end_ = end_points_[i];
      TRACE("Bounds are " << min_lens_[i] << " " << max_len_);
      TRACE("Current end vertex " << g_.int_id(current_end_));
      error_code |= Go(current_end_, min_lens_[i], 0, dijkstra);
      callback_.Flush();
    }
    return error_code; // 3 two mistakes, 2 bad dijkstra, 1 bad dfs, 0 = okay
  }

 private:
  static const size_t MAX_CALL_CNT          = 3000;
  static const size_t MAX_DIJKSTRA_VERTICES = 3000;

  const Graph& g_;
  const vector<size_t>& min_lens_;
  size_t max_len_;
  VertexId start_;
  const vector<VertexId>& end_points_;
  Callback& callback_;
  Path path_;
  size_t call_cnt_;

  DECL_LOGGER("PathProcessor")
};

template<class Graph>
class CompositeCallback : public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<EdgeId> Path;

 public:
  void AddProcessor(typename PathProcessor<Graph>::Callback & processor) {
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
class PathStorageCallback : public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<EdgeId> Path;

 public:
  PathStorageCallback(const Graph& g) : g_(g)
  {
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
class NonEmptyPathCounter : public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<EdgeId> Path;

 public:
  NonEmptyPathCounter(const Graph& g) : g_(g), count_(0)
  {
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
class VertexLabelerCallback : public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef vector<EdgeId> Path;

 public:
  VertexLabelerCallback(const Graph& g) : g_(g), count_(0)
  {
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
class DistancesLengthsCallback : public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<EdgeId> Path;

 public:
  DistancesLengthsCallback(const Graph& g) : g_(g)
  {
  }

  virtual void Flush() {
    all_distances_.push_back(distances_);          
    distances_.clear();
  }

  virtual void HandleReversedPath(const Path& path) {
    size_t path_length = 0;
    for (auto it = path.begin(); it != path.end(); ++it) {
      path_length += g_.length(*it);
    }
    distances_.push_back(path_length);
  }

  vector<size_t> distances(size_t k = 0) const {
    return all_distances_[k];
  }

 private:
  size_t PathLength(const Path& path) const {
    size_t res = 0;
    for (auto I = path.begin(); I != path.end(); ++I)
      res += g_.length(*I);
    return res;
  }

  const Graph& g_;
  vector<size_t> distances_;
  vector<vector<size_t>> all_distances_;
};

}

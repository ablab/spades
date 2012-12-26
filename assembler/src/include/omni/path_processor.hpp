#pragma once

#include "splitters.hpp"

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
public:
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;

  class Callback {
  public:
    virtual ~Callback() {

    }

    virtual void HandlePath(const vector<EdgeId>& path) = 0;
  };

private:
  const Graph& g_;
  size_t min_length_;
  size_t max_length_;
  VertexId start_;
  VertexId end_;
  Callback& callback_;

  vector<EdgeId> path_;

  size_t call_cnt_;

  static const size_t MAX_CALL_CNT = 30000;
  static const size_t MAX_DIJKSTRA_VERTICES = 4000;

//  void ReportLocality() {
//    TRACE("Backward dijkstra creation started");
//    BackwardBoundedDijkstra<Graph> backward_dijkstra(g_, max_length_);
//    TRACE("Backward dijkstra created with bound " << max_length_);
//    TRACE("Backward dijkstra started");
//    perf_counter pc;
//    backward_dijkstra.run(end_);
//    set<VertexId> visited;
//    stack<VertexId> to_go;
//    to_go.push(start_);
//    while (!visited.empty()) {
//      VertexId v = to_go.top();
//      to_go.pop();
//      visited.insert(v);
//      vector<EdgeId> out = g_.OutgoingEdges(v);
//
//      for (auto it = out.begin(); it != out.end(); ++it) {
//        VertexId end = g_.EdgeEnd(*it);
//        if (visited.count(end) == 0 && backward_dijkstra.DistanceCounted(end)) {
//          to_go.push(end);
//        }
//      }
//    }
//    ReliableSplitter<Graph/*, set<VertexId>::const_iterator*/> splitter(g_, visited.begin(), visited.end(), 50, 20000);
//
//    make_dir("br_strange_loc");
//    string folder = "br_strange_loc/" + ToString(g_.int_id(start_)) + "_" + ToString(g_.int_id(end_)) + "/";
//    make_dir(folder);
//
//    WriteComponents<Graph>(g_,
//        splitter,
//        folder + "graph.dot", *DefaultColorer(g_),
//        *StrGraphLabelerInstance(g_));
//  }

  //todo rewrite without recursion
  int Go(VertexId v, size_t current_path_length, Dijkstra<Graph>& distances_to_end) {
    int result = 0;
    TRACE("Processing vertex " << g_.int_id(v) << " started; current path length " << current_path_length);
    call_cnt_++;
    if (call_cnt_ == MAX_CALL_CNT) {
      DEBUG("Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
//      ReportLocality();
    }
    if (call_cnt_ >= MAX_CALL_CNT)
      return 1;

    if (!distances_to_end.DistanceCounted(v)
      || distances_to_end.GetDistance(v) + current_path_length > max_length_)
    {
      if (!distances_to_end.DistanceCounted(v)) {
        TRACE("Shortest distance from this vertex wasn't counted");
      } else if (distances_to_end.GetDistance(v) + current_path_length > max_length_) {
        TRACE("Shortest distance from this vertex is " << distances_to_end.GetDistance(v)
           << " and sum with current path length " << current_path_length
           << " exceeded max length " << max_length_);
      }
      return 0;
    }
    TRACE("Vertex " << g_.int_id(v) << " should be processed");

    if (v == end_ && current_path_length >= min_length_) {
      TRACE("New path found: " << PrintPath(g_, path_));
      callback_.HandlePath(path_);
      TRACE("Callback finished");
    }
    TRACE("Iterating through outgoing edges of vertex " << g_.int_id(v))
    for (auto I = g_.out_begin(v), E = g_.out_end(v); I != E; ++I) {
      EdgeId edge = *I;
      TRACE("Processing outgoing edge " << g_.int_id(edge) << " started");
      path_.push_back(edge);
      result |= Go(g_.EdgeEnd(edge), current_path_length + g_.length(edge), distances_to_end);
      path_.pop_back();
      TRACE("Processing outgoing edge " << g_.int_id(edge) << " finished");
    }
    TRACE("Processing vertex " << g_.int_id(v) << " finished");
    return result;
  }

public:
  PathProcessor(const Graph& g, double min_length, double max_length,
      VertexId start, VertexId end, Callback& callback) :
      g_(g), min_length_(
          (min_length < 0) ? 0 : (size_t) std::floor(min_length)), max_length_(
          (size_t) std::floor(max_length + 0.5)), start_(start), end_(
          end), callback_(callback), call_cnt_(0) {
    TRACE("Finding path from vertex " << g.int_id(start_) << " to vertex " << g.int_id(end_) << " of length [" << min_length_ << ", " << max_length_ << "]");
  }

  ~PathProcessor() {
    //WARN("Stopped looking for the path");
  }

  int Process() {
    int error_code = 0;
    TRACE("Backward dijkstra creation started");
    BackwardReliableBoundedDijkstra<Graph>
      backward_dijkstra(g_, max_length_, MAX_DIJKSTRA_VERTICES);
    TRACE("Backward dijkstra created with bound " << max_length_);

    backward_dijkstra.run(end_);

    if (backward_dijkstra.VertexLimitExceeded()) {
      TRACE("backward_dijkstra : vertex limit exceeded");
      error_code = 2;
    }

    TRACE("Starting recursive traversal");
    error_code |= Go(start_, 0, backward_dijkstra);
    if (call_cnt_ > 10)
      TRACE("number of calls: " << call_cnt_);
    TRACE("Recursive traversal finished");
    return error_code; // 3 two mistakes, 2 bad dijkstra, 1 bad dfs, 0 = okay
  }

private:
  DECL_LOGGER("PathProcessor")
};

template<class Graph>
class CompositeCallback: public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;

  vector<typename PathProcessor<Graph>::Callback*> processors_;
public:
  CompositeCallback(/*vector<PathProcessor<Graph>*>*/) {

  }

  virtual ~CompositeCallback() {

  }

  void AddProcessor(typename PathProcessor<Graph>::Callback & processor) {
    processors_.push_back(&processor);
  }

  virtual void HandlePath(const vector<EdgeId>& path) {
    for (auto it = processors_.begin(); it != processors_.end(); ++it) {
      (*it)->HandlePath(path);
    }
  }

};

template<class Graph>
class PathStorageCallback: public PathProcessor<Graph>::Callback {
public:
  typedef typename Graph::EdgeId EdgeId;

private:
  const Graph& g_;

  std::vector<vector<EdgeId>> paths_;
public:

  PathStorageCallback(const Graph& g) :
      g_(g) {
  }

  virtual void HandlePath(const vector<EdgeId>& path) {
    paths_.push_back(path);
  }

  size_t size() {
    return paths_.size();
  }

  std::vector<vector<EdgeId>> paths() {
    return paths_;
  }
};

template<class Graph>
class NonEmptyPathCounter: public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;

  //todo temporary
  const Graph& g_;

  size_t count_;
  vector<vector<EdgeId> > paths_;
public:

  NonEmptyPathCounter(const Graph& g) :
      g_(g), count_(0) {
    //    cout << "........................" << endl;
  }

  virtual void HandlePath(const vector<EdgeId>& path) {
    if (path.size() > 0) {
      //WARN("here " << path);

      /*
       size_t s = 0;
       for (auto it = path.begin(); it != path.end(); ++it) {
       s += g_.length(*it);
       cout << *it << "(" << g_.length(*it) << "), ";
       }

       cout << "Length " << s << endl;
       cout << "\n" << endl;
       */

      count_++;
      paths_.push_back(path);
    }
  }

  size_t count() {
    return count_;
  }
  vector<vector<EdgeId> > paths() {
    return paths_;
  }
};

template<class Graph>
class VertexLablerCallback: public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;

  Graph& g_;
  size_t count_;
  set<VertexId> vertices_;
public:

  VertexLablerCallback(Graph& g) :
      g_(g), count_(0) {
  }

  virtual void HandlePath(const vector<EdgeId>& path) {
    for (auto it = path.begin(); it != path.end(); ++it) {
      if (path.size() > 0) {
        vertices_.insert(g_.EdgeStart(*it));
        vertices_.insert(g_.EdgeEnd(*it));
        count_++;
      }
    }
  }

  const set<VertexId>& vertices() {
    return vertices_;
  }

  size_t count() {
    return count_;
  }
};

template<class Graph>
class DifferentDistancesCallback: public PathProcessor<Graph>::Callback {
  typedef typename Graph::EdgeId EdgeId;

  const Graph& g_;
  set<size_t> distances_;

public:
  DifferentDistancesCallback(const Graph& g) : g_(g)
  {
  }

  virtual ~DifferentDistancesCallback()
  {
  }

  virtual void HandlePath(const vector<EdgeId>& path) {
    size_t path_length = 0;
    for (auto it = path.begin(); it != path.end(); ++it) {
      path_length += g_.length(*it);
    }
    distances_.insert(path_length);
  }

  vector<size_t> distances() {
    return vector<size_t>(distances_.begin(), distances_.end());
  }

};
}

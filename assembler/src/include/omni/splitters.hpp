//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "dijkstra.hpp"
//#include "observable_graph.hpp"
namespace omnigraph {

template<class Element>
class ComponentSplitter {
public:
  virtual vector<Element> NextComponent() = 0;

  virtual bool Finished() = 0;

  virtual string ComponentName() {
    return "";
  }

  virtual ~ComponentSplitter() {
  }
};

template<class Element>
class PrecountedComponentSplitter : public ComponentSplitter<Element> {
  bool finished_;
  vector<Element> elements_;
public:

  template<class It>
  PrecountedComponentSplitter(It begin, It end) :
      finished_(false), elements_(begin, end) {

  }

  virtual vector<Element> NextComponent() {
    finished_ = true;
    return elements_;
  }

  virtual bool Finished() {
    return finished_;
  }

  virtual string ComponentName() {
    return "";
  }

};

template<class Graph>
class GraphSplitter: public ComponentSplitter<typename Graph::VertexId> {
  const Graph& graph_;
protected:

  GraphSplitter(const Graph& graph) :
      graph_(graph) {

  }

  const Graph& graph() const {
    return graph_;
  }
};

template<class Graph>
class ComponentFinder: public UnorientedDijkstra<Graph, size_t> {
private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef UnorientedDijkstra<Graph, size_t> base;
  set<EdgeId>& edges_;

public:
  ComponentFinder(const Graph &g, set<EdgeId> &edges) :
      base(g), edges_(edges) {
  }

  virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, size_t length) const {
    return edges_.count(edge) != 0;
  }
};

template<class Graph>
class NeighbourhoodFinder: public UnorientedDijkstra<Graph, size_t> {
private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef UnorientedDijkstra<Graph, size_t> base;
  set<EdgeId>& edges_;
  const size_t bound_;

public:
  NeighbourhoodFinder(const Graph &g, set<EdgeId> &edges, size_t bound) :
      base(g), edges_(edges), bound_(bound) {
  }

  virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
    return distance <= bound_;
  }

  virtual size_t GetLength(EdgeId edge) const {
    if (edges_.count(edge) != 0)
      return 0;
    else
      return this->graph().length(edge);
  }

};

template<class Graph>
class SubgraphDijkstra: public UnorientedDijkstra<Graph, size_t> {
private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef UnorientedDijkstra<Graph, size_t> base;
  const set<VertexId> &subgraph_;

public:
  SubgraphDijkstra(const Graph& g, const set<VertexId>& subgraph) :
      base(g), subgraph_(subgraph) {
  }

  virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, size_t length) const {
    return subgraph_.count(vertex) != 0;
  }

};

template<class Graph>
class ErrorComponentSplitter: public GraphSplitter<Graph> {
private:
  typedef GraphSplitter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  set<EdgeId> black_edges_;
  typename Graph::SmartEdgeIt iterator_;
  set<VertexId> visited_;

public:
  ErrorComponentSplitter(const Graph &graph, const set<EdgeId> &black_edges) :
      base(graph), black_edges_(black_edges), iterator_(
          graph.SmartEdgeBegin()) {
    TRACE("ErrorComponentSplitter created and SmartIterator initialized");
  }

  virtual ~ErrorComponentSplitter() {
  }

  vector<VertexId> FindComponent(VertexId start_vertex) {
    ComponentFinder<Graph> cf(this->graph(), black_edges_);
    cf.run(start_vertex);
    return cf.ReachedVertices();
  }

  vector<VertexId> FindNeighbourhood(VertexId start, size_t bound) {
    NeighbourhoodFinder<Graph> nf(this->graph(), black_edges_, bound);
    nf.run(start);
    return nf.ReachedVertices();
  }

  size_t FindDiameter(const vector<VertexId> &component) {
    set < VertexId > component_set(component.begin(), component.end());
    size_t result = 0;
    VertexId current = *(component.begin());
    for (size_t i = 0; i < 4; i++) {
      pair<VertexId, size_t> next = GetFarthest(current, component_set);
      current = next.first;
      result = next.second;
    }
    return result;
  }

  pair<VertexId, size_t> GetFarthest(VertexId v,
      const set<VertexId> &component) {
    SubgraphDijkstra<Graph> sd(this->graph(), component);
    sd.run(v);
    pair<VertexId, size_t> result(v, 0);
    auto bounds = sd.GetDistances();
    for (auto it = bounds.first; it != bounds.second; ++it) {
      if (it->second > result.second) {
        result = *it;
      }
    }
    return result;
  }

  virtual vector<VertexId> NextComponent() {
    TRACE("Construction of next component started");
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    EdgeId next = *iterator_;
    ++iterator_;
    vector < VertexId > component = FindComponent(
        this->graph().EdgeEnd(next));
    TRACE(
        "Error edges component constructed. It contains "
            << component.size() << " vertices");
    size_t component_size = FindDiameter(component);
    TRACE("Diameter of component is " << component_size);
    vector < VertexId > neighbourhood = FindNeighbourhood(
        this->graph().EdgeEnd(next), 1.5 * component_size);
    TRACE(
        "Error edges component neighborhood constructed. It contains "
            << neighbourhood.size() << " vertices");
    visited_.insert(component.begin(), component.end());
    return neighbourhood;
  }

  virtual bool Finished() {
    while (!iterator_.IsEnd()) {
      if (black_edges_.find(*iterator_) != black_edges_.end()
          && visited_.find(this->graph().EdgeEnd(*iterator_))
              == visited_.end()) {
        return false;
      }
      ++iterator_;
    }
    return true;
  }

};

template<class Graph>
class ShortEdgeComponentNeighbourhoodFinder: public UnorientedDijkstra<Graph> {
private:
  typedef UnorientedDijkstra<Graph> base;
protected:
  typedef typename base::VertexId VertexId;
  typedef typename base::EdgeId EdgeId;
  typedef typename base::DistanceType distance_t;
private:
  distance_t bound_;
public:
  ShortEdgeComponentNeighbourhoodFinder(const Graph &graph, distance_t bound) :
      UnorientedDijkstra<Graph>(graph), bound_(bound) {
  }

  virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
    return distance == 0;
  }

  virtual distance_t GetLength(EdgeId edge) const {
    if (this->graph().length(edge) <= bound_)
      return 0;
    else
      return 1;
  }
};

template<class Graph>
class LongEdgesInclusiveSplitter: public GraphSplitter<Graph> {
private:
  typedef GraphSplitter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  erasable_priority_queue<VertexId, std::less<VertexId>> queue_;
  //	SmartVertexIterator<omnigraph::ObservableGraph<VertexId, EdgeId> >
  //			iterator_;
  set<VertexId> visited_;
  size_t bound_;

public:
  LongEdgesInclusiveSplitter(const Graph &graph, size_t bound) :
      base(graph), queue_(graph.begin(), graph.end()), /*iterator_(graph.SmartVertexBegin()), */
      bound_(bound) {
  }

  virtual vector<VertexId> NextComponent() {
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    VertexId next = queue_.top();
    TRACE("Search started");
    queue_.pop();
    ShortEdgeComponentNeighbourhoodFinder<Graph> cf(this->graph(), bound_);
    cf.run(next);
    TRACE("Search finished");
    vector < VertexId > result = cf.ReachedVertices();
    for (auto it = result.begin(); it != result.end(); ++it) {
      if (cf.GetDistance(*it) == 0) {
        //				iterator_.erase(*it);
        queue_.erase(*it);
      }
    }
    TRACE("Component vector filled");
    return result;
  }

  virtual bool Finished() {
    //		return iterator_.IsEnd();
    return queue_.empty();
  }

};

template<class Graph, typename distance_t = size_t>
class CountingDijkstra: public UnorientedDijkstra<Graph, distance_t> {
private:
  typedef UnorientedDijkstra<Graph, distance_t> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  static const distance_t inf = 100000000;

  const size_t max_size_;

  const size_t edge_length_bound_;

  mutable size_t current_;

public:
  CountingDijkstra(const Graph &graph, size_t max_size,
      size_t edge_length_bound) :
      base(graph), max_size_(max_size), edge_length_bound_(
          edge_length_bound), current_(0) {
  }

  virtual bool CheckPutVertex(VertexId vertex, EdgeId edge,
      distance_t length) const {
    if (current_ < max_size_) {
      ++current_;
    }
    if (current_ < max_size_ && GetLength(edge) < inf) {
      return true;
    }
    return false;
  }

  virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
    return current_ < max_size_;
  }

  virtual void init(VertexId start) {
    current_ = 0;
  }

  virtual size_t GetLength(EdgeId edge) const {
    if (this->graph().length(edge) <= edge_length_bound_)
      //todo change back
//			return 1;
      return this->graph().length(edge);
    else
      return inf;
  }
};

template<class Graph, typename distance_t = size_t>
class CountingDijkstraForPaths: public CountingDijkstra<Graph, distance_t> {
  typedef CountingDijkstra<Graph, distance_t> base;
protected:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
private:
  set<VertexId> path_vertices_;

  set<VertexId> CollectVertices(const vector<EdgeId>& path) {
    set < VertexId > answer;
    for (auto it = path.begin(); it != path.end(); ++it) {
      answer.insert(this->graph().EdgeStart(*it));
      answer.insert(this->graph().EdgeEnd(*it));
    }
    return answer;
  }

public:
  CountingDijkstraForPaths(const Graph& graph, size_t max_size,
      size_t edge_length_bound, const vector<EdgeId>& path) :
      base(graph, max_size, edge_length_bound), path_vertices_(
          CollectVertices(path)) {
  }

  virtual size_t GetLength(EdgeId edge) const {
    if (path_vertices_.count(this->graph().EdgeStart(edge))
        && path_vertices_.count(this->graph().EdgeEnd(edge)))
        return min(int(base::GetLength(edge)), 200);
//			return 1;
    return base::GetLength(edge);
//		return 2;
  }

};

template<class Graph, typename distance_t = size_t>
class ComponentCloser {
private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  const Graph &graph_;
  size_t edge_length_bound_;

public:
  ComponentCloser(const Graph &graph, size_t edge_length_bound) :
      graph_(graph), edge_length_bound_(edge_length_bound) {
  }

  ~ComponentCloser() {
  }

  void AddNewVertices(vector<EdgeId> possible_long_edges,
      set<VertexId> &additional_vertices) {
  }

  void CloseComponent(set<VertexId> &component) {
    set<VertexId> additional_vertices;
    for (auto it = component.begin(); it != component.end(); ++it) {
      {
        vector<EdgeId> possible_long_edges = graph_.OutgoingEdges(*it);
        for (auto it = possible_long_edges.begin();
            it != possible_long_edges.end(); ++it) {
          if (graph_.length(*it) >= edge_length_bound_) {
            additional_vertices.insert(graph_.EdgeEnd(*it));
          }
        }
      }
      {
        vector<EdgeId> possible_long_edges = graph_.IncomingEdges(*it);
        for (auto it = possible_long_edges.begin();
            it != possible_long_edges.end(); ++it) {
          if (graph_.length(*it) >= edge_length_bound_) {
            additional_vertices.insert(graph_.EdgeStart(*it));
          }
        }
      }
    }
    component.insert(additional_vertices.begin(),
        additional_vertices.end());
  }
};

template<class Graph>
class ShortEdgeComponentFinder: public UnorientedDijkstra<Graph> {
  typedef UnorientedDijkstra<Graph> base;
  typedef typename base::DistanceType distance_t;
  typedef typename base::VertexId VertexId;
  typedef typename base::EdgeId EdgeId;

  distance_t bound_;
public:
  ShortEdgeComponentFinder(const Graph &graph, distance_t bound) :
      base(graph), bound_(bound) {
  }

  virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
    return distance == 0;
  }

  virtual distance_t GetLength(EdgeId edge) const {
    if (this->graph().length(edge) <= bound_)
      return 0;
    else
      return 1;
  }
};

template<class Graph, class It = typename Graph::VertexIterator>
class ReliableSplitter: public GraphSplitter<Graph> {
private:
  typedef GraphSplitter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  set<VertexId> visited_;
  size_t max_size_;
  size_t edge_length_bound_;
  It current_;
  It end_;

  void SkipVisited() {
    while (current_ != end_ && visited_.count(*current_) == 1) {
      ++current_;
    }
  }

public:
  ReliableSplitter(const Graph &graph, It begin, It end, size_t max_size,
      size_t edge_length_bound) :
      base(graph), max_size_(max_size), edge_length_bound_(
          edge_length_bound), current_(begin), end_(end) {
    TRACE(
        "Long edges splitter created and queue filled with all graph vertices");
  }

  ReliableSplitter(const Graph &graph, size_t max_size,
      size_t edge_length_bound) :
      base(graph), max_size_(max_size), edge_length_bound_(
          edge_length_bound), current_(graph.begin()), end_(
          graph.end()) {
    TRACE(
        "Long edges splitter created and queue filled with all graph vertices");
  }

  virtual ~ReliableSplitter() {
  }

  virtual vector<VertexId> NextComponent() {
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    TRACE("Search started");
    CountingDijkstra<Graph> cf(this->graph(), max_size_,
        edge_length_bound_);
    cf.run(*current_);
    TRACE("Search finished");
    //TODO Refactor this!!!!
    vector < VertexId > result_vector = cf.ReachedVertices();
    set < VertexId > result(result_vector.begin(), result_vector.end());
//		set<VertexId, typename Graph::Comparator> result(cf.ProcessedVertices());
    visited_.insert(result.begin(), result.end());
    ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
    cc.CloseComponent(result);
    TRACE("Component vector filled");
    SkipVisited();
    return vector<VertexId>(result.begin(), result.end());
  }

  virtual bool Finished() {
    return current_ == end_;
  }
};

template<class Graph>
class ReliableSplitterAlongPath: public GraphSplitter<Graph> {
private:
  typedef GraphSplitter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  size_t max_size_;
  size_t edge_length_bound_;
  set<VertexId> last_component_;
  size_t current_index_;
  MappingPath<EdgeId> path_;
  Range covered_range_;
  bool start_processed_;

  //todo edge not used in the body
  bool EdgeCovered(EdgeId edge) {
    return last_component_.count(
        this->graph().EdgeStart(path_[current_index_].first)) == 1
        && last_component_.count(
            this->graph().EdgeEnd(path_[current_index_].first)) == 1;
  }

  void SkipVisited() {
    covered_range_.start_pos =
        path_[current_index_].second.initial_range.start_pos;
    covered_range_.end_pos =
        path_[current_index_].second.initial_range.end_pos;
    //always go forward at least one path element
    ++current_index_;
    while (current_index_ < path_.size()
        && EdgeCovered(path_[current_index_].first)) {
      covered_range_.end_pos =
          path_[current_index_].second.initial_range.end_pos;
      ++current_index_;
    }
  }

public:
  string ComponentName() {
    stringstream ss;
    ss << covered_range_.start_pos << "_" << covered_range_.end_pos;
    return ss.str();
  }

  ReliableSplitterAlongPath(const Graph &graph, size_t max_size,
      size_t edge_length_bound, const MappingPath<EdgeId>& path) :
      base(graph), max_size_(max_size), edge_length_bound_(
          edge_length_bound), current_index_(0), path_(path), covered_range_(
          0, 0), start_processed_(false) {

  }

  virtual ~ReliableSplitterAlongPath() {
  }

  virtual vector<VertexId> NextComponent() {
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    TRACE("Search started");
    CountingDijkstraForPaths<Graph> cf(this->graph(), max_size_,
        edge_length_bound_, path_.simple_path().sequence());
    if (start_processed_)
      cf.run(this->graph().EdgeEnd(path_[current_index_].first));
    else {
      cf.run(this->graph().EdgeStart(path_[current_index_].first));
      start_processed_ = true;
    }
    TRACE("Search finished");
    set < VertexId > last_component = cf.ProcessedVertices();
    last_component_.clear();
    last_component_.insert(
        this->graph().EdgeStart(path_[current_index_].first));
    last_component_.insert(
        this->graph().EdgeEnd(path_[current_index_].first));
    last_component_.insert(last_component.begin(), last_component.end());

    TRACE("Component vector filled");
    size_t prev_index = current_index_;
    ComponentCloser<Graph> cc(this->graph(), 1/*edge_length_bound_*/);
    cc.CloseComponent(last_component_);
    SkipVisited();
    //todo ask Anton what is this...
    if (prev_index + 1 != current_index_) {
      start_processed_ = true;
    } else if (!start_processed_) {
      current_index_ = prev_index;
      start_processed_ = true;
    } else {
      start_processed_ = false;
    }
    return vector<VertexId>(last_component_.begin(), last_component_.end());
  }

  virtual bool Finished() {
    return current_index_ == path_.size();
  }
};

template<class Graph>
class LongEdgesExclusiveSplitter: public GraphSplitter<Graph> {
private:
  typedef GraphSplitter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
//	erasable_priority_queue<VertexId> queue_;
  typename Graph::SmartVertexIt iterator_;
  set<VertexId> visited_;
  size_t bound_;
private:
  DECL_LOGGER("LongEdgesExclusiveSplitter")
public:
  LongEdgesExclusiveSplitter(const Graph &graph, size_t bound) :
  base(graph), /*queue_(graph.begin(), graph.end()), */iterator_(graph.SmartVertexBegin()),
  bound_(bound) {
  }

  virtual vector<VertexId> NextComponent() {
    TRACE("search started");
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    VertexId next = *iterator_; //queue_.top();
    ++iterator_; //queue_.pop();
    ShortEdgeComponentFinder<Graph> cf(this->graph(), bound_);
    cf.run(next);

    TRACE("comp Finder finished");
    auto result = cf.ProcessedVertices();
    for (auto it = result.begin(); it != result.end(); ++it) {
      iterator_.erase(*it); //queue_.erase(*it);
    }
    TRACE("Returning component");
    return vector<VertexId>(result.begin(), result.end());
  }

  virtual bool Finished() {
    TRACE("Inside Finished");
    return iterator_.IsEnd();
//		return queue_.empty();
  }

};

template<class Element>
class AbstractFilter {
public:
  virtual ~AbstractFilter() {
  }

  virtual bool Check(const Element &element) const = 0;
};

template<class Element>
class TrueFilter: public AbstractFilter<Element> {

  /*virtual*/
  bool Check(const Element &element) const {
    return true;
  }
};

template<class Graph>
class GraphComponentFilter: public AbstractFilter<
    vector<typename Graph::VertexId>> {
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  const Graph& graph_;
protected:
  GraphComponentFilter(const Graph& graph) :
      graph_(graph) {

  }

  const Graph& graph() const {
    return graph_;
  }
};

template<class Graph>
class AnyEdgeContainFilter: public GraphComponentFilter<Graph> {
  typedef GraphComponentFilter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  const vector<EdgeId> edges_of_interest_;
public:
  AnyEdgeContainFilter(const Graph& graph,
      const vector<EdgeId>& edges_of_interest) :
      base(graph), edges_of_interest_(edges_of_interest) {

  }

  AnyEdgeContainFilter(const Graph& graph, EdgeId edge_of_interest) :
      base(graph), edges_of_interest_( { edge_of_interest }) {

  }

  bool ContainsEdge(const vector<VertexId>& component, EdgeId e) const {
    return std::find(component.begin(), component.end(),
        this->graph().EdgeStart(e)) != component.end()
        && std::find(component.begin(), component.end(),
            this->graph().EdgeEnd(e)) != component.end();
  }

  /*virtual*/
  bool Check(const vector<VertexId> &component) const {
    for (auto it = edges_of_interest_.begin();
        it != edges_of_interest_.end(); ++it) {
      if (ContainsEdge(component, *it)) {
        return true;
      }
    }
    return false;
  }

};

template<class Graph>
class ComponentSizeFilter: public GraphComponentFilter<Graph> {
private:
  typedef GraphComponentFilter<Graph> base;
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  size_t max_length_;
  size_t min_vertex_number_;
  size_t max_vertex_number_;
public:
  ComponentSizeFilter(const Graph &graph, size_t max_length,
      size_t min_vertex_number, size_t max_vertex_number) :
      base(graph), max_length_(max_length), min_vertex_number_(min_vertex_number),
      max_vertex_number_(max_vertex_number) {
  }

  /*virtual*/
  bool Check(const vector<VertexId> &vertices) const {
    if (vertices.size() <= min_vertex_number_ || vertices.size() >= max_vertex_number_)
      return false;
    set < VertexId > component(vertices.begin(), vertices.end());
    for (auto iterator = vertices.begin(); iterator != vertices.end();
        ++iterator) {
      vector < EdgeId > edges = this->graph().OutgoingEdges(*iterator);
      for (auto edge_iterator = edges.begin();
          edge_iterator != edges.end(); edge_iterator++) {
        if (component.count(this->graph().EdgeEnd(*edge_iterator)) == 1
            && this->graph().length(*edge_iterator)
                <= max_length_) {
          return true;
        }
      }
    }
    return false;
  }
};

template<class Graph>
class FilteringSplitterWrapper: public ComponentSplitter<
    typename Graph::VertexId> {
private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

  ComponentSplitter<typename Graph::VertexId> &inner_splitter_;
  vector<VertexId> next;
  const AbstractFilter<vector<VertexId>> &checker_;
  bool ready;
public:
  FilteringSplitterWrapper(
      ComponentSplitter<typename Graph::VertexId> &inner_splitter,
      const AbstractFilter<vector<VertexId>> &checker) :
      inner_splitter_(inner_splitter), checker_(checker), ready(false) {
  }

  /*virtual*/
  string ComponentName() {
    return inner_splitter_.ComponentName();
  }

  /*virtual*/
  vector<VertexId> NextComponent() {
    if (Finished()) {
      VERIFY(false);
      return vector<VertexId>();
    }
    ready = false;
    return next;
  }

  /*virtual*/
  bool Finished() {
    if (!ready) {
      TRACE("Calculating next nontrivial component");
      while (!inner_splitter_.Finished()) {
        TRACE("Calculating next component");
        next = inner_splitter_.NextComponent();
        TRACE("Next component calculated");
        if (checker_.Check(next)) {
          TRACE("Nontrivial component found");
          ready = true;
          return false;
        }
        TRACE("Component skipped");
      }
      return true;
    }
    return false;
  }

private:
  DECL_LOGGER("FilteringSplitterWrapper");
};

template<class Graph>
class VertexNeighborhoodFinder: public GraphSplitter<Graph> {
private:
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  VertexId vertex_;
  size_t max_size_;
  size_t edge_length_bound_;
  bool finished_;
public:
  VertexNeighborhoodFinder(const Graph &graph, VertexId vertex,
      size_t max_size, size_t edge_length_bound) :
      GraphSplitter<Graph>(graph), vertex_(vertex), max_size_(max_size), edge_length_bound_(
          edge_length_bound), finished_(false) {
  }

  /*virtual*/
  vector<VertexId> NextComponent() {
    CountingDijkstra<Graph> cf(this->graph(), max_size_,
        edge_length_bound_);
    set < VertexId > result_set;
    cf.run(vertex_);
    vector < VertexId > result = cf.ReachedVertices();
    result_set.insert(result.begin(), result.end());

    ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
    cc.CloseComponent(result_set);

    finished_ = true;
    return vector<VertexId>(result_set.begin(), result_set.end());
  }

  /*virtual*/
  bool Finished() {
    return finished_;
  }
};

template<class Graph>
class EdgeNeighborhoodFinder: public GraphSplitter<Graph> {
private:
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  EdgeId edge_;
  size_t max_size_;
  size_t edge_length_bound_;
  bool finished_;
public:
  EdgeNeighborhoodFinder(const Graph &graph, EdgeId edge, size_t max_size,
      size_t edge_length_bound) :
      GraphSplitter<Graph>(graph), edge_(edge), max_size_(max_size), edge_length_bound_(
          edge_length_bound), finished_(false) {
  }

  /*virtual*/
  vector<VertexId> NextComponent() {
    CountingDijkstra<Graph> cf(this->graph(), max_size_,
        edge_length_bound_);
    set < VertexId > result_set;
    cf.run(this->graph().EdgeStart(edge_));
    vector < VertexId > result_start = cf.ReachedVertices();
    result_set.insert(result_start.begin(), result_start.end());
    cf.run(this->graph().EdgeEnd(edge_));
    vector < VertexId > result_end = cf.ReachedVertices();
    result_set.insert(result_end.begin(), result_end.end());

    ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
    cc.CloseComponent(result_set);

    finished_ = true;
    vector < VertexId > result;
    for (auto it = result_set.begin(); it != result_set.end(); ++it)
      result.push_back(*it);
    return result;
  }

  /*virtual*/
  bool Finished() {
    return finished_;
  }
};
}

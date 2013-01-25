//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef DIJKSTRA_HPP_
#define DIJKSTRA_HPP_

#include <queue>
#include <map>
#include <set>

namespace omnigraph {

template<typename distance_t, typename T>
class ReverseDistanceComparator {
public:
  ReverseDistanceComparator() {
  }

  bool operator()(pair<distance_t, T> a, pair<distance_t, T> b) {
    if(a.first != b.first)
      return b.first < a.first;
    else
      return b.second < a.second;
  }
};

template<class Graph, typename distance_t = size_t>
class Dijkstra {
protected:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef distance_t DistanceType;
  typedef std::map<VertexId, distance_t> distances_map;
  typedef typename distances_map::const_iterator distances_map_ci;
  typedef std::pair<distance_t, VertexId> element_t;
  typedef typename std::priority_queue<element_t, vector<element_t>,
                                       ReverseDistanceComparator<distance_t, VertexId>> queue_t;

 public:
  bool finished() const {
    return finished_;
  }

  Dijkstra(const Graph& graph) :
      graph_(graph), finished_(false)
  {
  }

  virtual ~Dijkstra()
  {
  }

  bool DistanceCounted(VertexId vertex) const {
    return distances_.find(vertex) != distances_.end();
  }

  virtual void init(VertexId start) {
  }

  distance_t GetDistance(VertexId vertex) const {
    VERIFY(DistanceCounted(vertex));
    return distances_.find(vertex)->second;
  }

  pair<distances_map_ci, distances_map_ci> GetDistances() const {
    distances_map_ci begin = distances_.begin();
    distances_map_ci end = distances_.end();
    return make_pair(begin, end);
  }

  void set_finished(bool state) {
    finished_ = state;
  }

  virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) const {
    return true;
  }

  virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
    return true;
  }

  virtual distance_t GetLength(EdgeId edge) const {
    return graph_.length(edge);
  }

  virtual void AddNeighboursToQueue(VertexId cur_vertex, distance_t cur_dist, queue_t& queue) {
    BOOST_FOREACH(EdgeId edge, graph_.OutgoingEdges(cur_vertex)) {
    //for (auto I = graph_.out_begin(cur_vertex), E = graph_.out_end(cur_vertex); I != E; ++I) {
      //EdgeId edge = *I;
      //TRACE("Checking " << i << "th neighbour of vertex " << graph_.str(cur_vertex) << " started");
      VertexId neighbour = graph_.EdgeEnd(edge);
      if (!DistanceCounted(neighbour)) {
        //TRACE("Adding new entry to queue");
        distance_t new_dist = GetLength(edge) + cur_dist;
        //TRACE("Entry: vertex " << graph_.str(cur_vertex) << " distance " << new_distance);
        if (CheckPutVertex(neighbour, edge, new_dist)) {
        //TRACE("CheckPutVertex returned true and new entry is added");
          queue.push(make_pair(new_dist, neighbour));
        }
      }
      //TRACE("Checking " << i << "th neighbour of vertex " << graph_.str(cur_vertex) << " finished");
    }
    TRACE("All neighbours of vertex " << graph_.str(cur_vertex) << " processed");
  }

  void run(VertexId start) {
    TRACE("Starting dijkstra run from vertex " << graph_.str(start));
    set_finished(false);
    distances_.clear();
    processed_vertices_.clear();
    init(start);
    queue_t queue;
    queue.push(make_pair(0, start));
    TRACE("Priority queue initialized. Starting search");

    while (!queue.empty() && !finished()) {
      TRACE("Dijkstra iteration started");
      const element_t& next = queue.top();
      distance_t distance = next.first;
      VertexId vertex = next.second;
      queue.pop();
      //TRACE("Vertex " << graph_.str(vertex) << " with distance " << distance << " fetched from queue");
      if (DistanceCounted(vertex)) {
        //TRACE("Distance to vertex " << graph_.str(vertex) << " already counted. Proceeding to next queue entry.");
        continue;
      }
      distances_.insert(make_pair(vertex, distance));
      //TRACE("Vertex " << graph_.str(vertex) << " is found to be at distance "
              //<< distance << " from vertex " << graph_.str(start));
      if (!CheckProcessVertex(vertex, distance)) {
        //TRACE("Check for processing vertex failed. Proceeding to the next queue entry.");
        continue;
      }

      processed_vertices_.insert(vertex);
      this->AddNeighboursToQueue(vertex, distance, queue);
    }
    set_finished(true);
    TRACE("Finished dijkstra run from vertex " << graph_.str(start));
  }

  vector<VertexId> ReachedVertices() const {
    vector<VertexId> result;
    for (auto it = distances_.begin(); it != distances_.end(); ++it) {
      result.push_back(it->first);
    }
    return result;
  }

  const set<VertexId>& ProcessedVertices() {
    return processed_vertices_;
  }

 protected:
  const Graph& graph() const {
    return graph_;
  }

 private:
  const Graph& graph_;
  bool finished_;
  distances_map distances_;
  set<VertexId> processed_vertices_;

  DECL_LOGGER("Dijkstra");
};

template<class Graph>
class DistanceCounter {
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;

public:
  DistanceCounter(const Graph& graph) :
    graph_(graph), dijkstra_(graph), ready_(false) {
  }

  bool IsReachable(VertexId from, VertexId to) {
    EnsureFrom(from);
    return dijkstra_.DistanceCounted(to);
  }

  size_t Distance(VertexId from, VertexId to) {
    EnsureFrom(from);
    return dijkstra_.GetDistance(to);
  }

private:
  void EnsureFrom(VertexId from) {
    if (!ready_ || prev_ != from) {
      dijkstra_.run(from);
      ready_ = true;
      prev_ = from;
    }
  }

  const Graph& graph_;
  Dijkstra<Graph> dijkstra_;
  VertexId prev_;
  bool ready_;
};

template<class Graph, typename distance_t = size_t>
class UnorientedDijkstra: public Dijkstra<Graph, distance_t> {
private:
  typedef Dijkstra<Graph, distance_t> base;
protected:
  typedef typename base::VertexId VertexId;
  typedef typename base::EdgeId EdgeId;
  typedef typename base::element_t element_t;
  typedef typename base::queue_t queue_t;

public:
  UnorientedDijkstra(const Graph &graph) : base(graph) {}

  virtual void AddNeighboursToQueue(VertexId cur_vertex, distance_t cur_dist, queue_t& queue) {
    const Graph& g = this->graph();
    for (auto I = g.out_begin(cur_vertex), E = g.out_end(cur_vertex); I != E; ++I) {
      EdgeId edge = *I;
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " started");
      VertexId neighbour = g.EdgeEnd(edge);
      if (!this->DistanceCounted(neighbour)) {
        //TRACE("Adding new entry to queue");
        distance_t new_dist = this->GetLength(edge) + cur_dist;
        //TRACE("Entry: vertex " << g.str(cur_vertex) << " distance " << new_distance);
        if (this->CheckPutVertex(neighbour, edge, new_dist)) {
        //TRACE("CheckPutVertex returned true and new entry is added");
          queue.push(make_pair(new_dist, neighbour));
        }
      }
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " finished");
    }

    for (auto I = g.in_begin(cur_vertex), E = g.in_end(cur_vertex); I != E; ++I) {
      EdgeId edge = *I;
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " started");
      VertexId neighbour = g.EdgeStart(edge);
      if (!this->DistanceCounted(neighbour)) {
        //TRACE("Adding new entry to queue");
        distance_t new_dist = this->GetLength(edge) + cur_dist;
        //TRACE("Entry: vertex " << g.str(cur_vertex) << " distance " << new_distance);
        if (this->CheckPutVertex(neighbour, edge, new_dist)) {
        //TRACE("CheckPutVertex returned true and new entry is added");
          queue.push(make_pair(new_dist, neighbour));
        }
      }
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " finished");
    }
    TRACE("All neighbours of vertex " << g.str(cur_vertex) << " processed");
  }
};

template<class Graph, typename distance_t = size_t>
class BackwardDijkstra: public Dijkstra<Graph, distance_t> {
  typedef Dijkstra<Graph, distance_t> base;
protected:
  typedef typename base::VertexId VertexId;
  typedef typename base::EdgeId EdgeId;
  typedef typename base::element_t element_t;
  typedef typename base::queue_t queue_t;

public:
  BackwardDijkstra(const Graph &graph) :
    base(graph) {
  }

  virtual void AddNeighboursToQueue(VertexId cur_vertex, distance_t cur_dist, queue_t& queue) {
    const Graph& g = this->graph();
    for (auto I = g.in_begin(cur_vertex), E = g.in_end(cur_vertex); I != E; ++I) {
      EdgeId edge = *I;
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " started");
      VertexId neighbour = g.EdgeStart(edge);
      if (!DistanceCounted(neighbour)) {
        //TRACE("Adding new entry to queue");
        distance_t new_dist = GetLength(edge) + cur_dist;
        //TRACE("Entry: vertex " << g.str(cur_vertex) << " distance " << new_distance);
        if (CheckPutVertex(neighbour, edge, new_dist)) {
        //TRACE("CheckPutVertex returned true and new entry is added");
          queue.push(make_pair(new_dist, neighbour));
        }
      }
      //TRACE("Checking " << i << "th neighbour of vertex " << g.str(cur_vertex) << " finished");
    }
    TRACE("All neighbours of vertex " << g.str(cur_vertex) << " processed");
  }

private:
  DECL_LOGGER("BackwardDijkstra");
};

template<class Graph, typename distance_t = size_t>
class BoundedDijkstra: public Dijkstra<Graph, distance_t> {
  typedef Dijkstra<Graph, distance_t> base;
protected:
  typedef typename base::VertexId VertexId;
  typedef typename base::EdgeId EdgeId;

public:
  BoundedDijkstra(const Graph &graph, distance_t bound) :
      base(graph), bound_(bound) {}

  virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) const {
    return (length <= bound_);
  }

  virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
    return (distance <= bound_);
  }

private:
  distance_t bound_;
};

}
#endif /* DIJKSTRA_HPP_ */

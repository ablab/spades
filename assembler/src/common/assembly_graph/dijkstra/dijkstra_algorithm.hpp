//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "utils/simple_tools.hpp"
#include "dijkstra_settings.hpp"

#include <queue>
#include <vector>
#include <set>
#include <map>

namespace omnigraph {

template<typename Graph, typename distance_t = size_t>
struct element_t{
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    distance_t distance;
    VertexId curr_vertex;
    VertexId prev_vertex;
    EdgeId edge_between;

    element_t(distance_t new_distance, VertexId new_cur_vertex, VertexId new_prev_vertex,
            EdgeId new_edge_between) : distance(new_distance), curr_vertex(new_cur_vertex),
                    prev_vertex(new_prev_vertex), edge_between(new_edge_between) { }
};

template<typename T>
class ReverseDistanceComparator {
public:
  ReverseDistanceComparator() {
  }

  bool operator()(T obj1, T obj2){
      if(obj1.distance != obj2.distance)
          return obj2.distance < obj1.distance;
      if(obj2.curr_vertex != obj1.curr_vertex)
          return obj2.curr_vertex < obj1.curr_vertex;
      if(obj2.prev_vertex != obj1.prev_vertex)
          return obj2.prev_vertex < obj1.prev_vertex;
      return obj2.edge_between < obj1.edge_between;
  }
};

template<class Graph, class DijkstraSettings, typename distance_t = size_t>
class Dijkstra {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef distance_t DistanceType;

    typedef std::map<VertexId, distance_t> distances_map;
    typedef typename distances_map::const_iterator distances_map_ci;
    typedef typename std::priority_queue<element_t<Graph, distance_t>, std::vector<element_t<Graph, distance_t>>,
            ReverseDistanceComparator<element_t<Graph, distance_t>>> queue_t;

    // constructor parameters
    const Graph& graph_;
    DijkstraSettings settings_;
    const size_t max_vertex_number_;

    // changeable parameters
    bool finished_;
    size_t vertex_number_;
    bool vertex_limit_exceeded_;

    // accumulative structures
    distances_map distances_;
    std::set<VertexId> processed_vertices_;
    std::map<VertexId, pair<VertexId, EdgeId>> prev_vert_map_;

    void Init(VertexId start, queue_t &queue) {
        vertex_number_ = 0;
        distances_.clear();
        processed_vertices_.clear();
        prev_vert_map_.clear();
        set_finished(false);
        settings_.Init(start);
        queue.push(element_t<Graph, distance_t>(0, start, VertexId(0), EdgeId(0)));
        prev_vert_map_[start] = std::pair<VertexId, EdgeId>(VertexId(0), EdgeId(0));
    }

    void set_finished(bool state) {
        finished_ = state;
    }

    bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) const {
        return settings_.CheckPutVertex(vertex, edge, length);
    }

    bool CheckProcessVertex(VertexId vertex, distance_t distance) {
        ++vertex_number_;
        if (vertex_number_ > max_vertex_number_) {
            vertex_limit_exceeded_ = true;
            return false;
        }
        return (vertex_number_ < max_vertex_number_) && settings_.CheckProcessVertex(vertex, distance);
    }

    distance_t GetLength(EdgeId edge) const {
        return settings_.GetLength(edge);
    }

    void AddNeighboursToQueue(VertexId cur_vertex, distance_t cur_dist, queue_t& queue) {
        auto neigh_iterator = settings_.GetIterator(cur_vertex);
        while (neigh_iterator.HasNext()) {
            TRACE("Checking new neighbour of vertex " << graph_.str(cur_vertex) << " started");
            auto cur_pair = neigh_iterator.Next();
            if (!DistanceCounted(cur_pair.vertex)) {
                TRACE("Adding new entry to queue");
                distance_t new_dist = GetLength(cur_pair.edge) + cur_dist;
                TRACE("Entry: vertex " << graph_.str(cur_vertex) << " distance " << new_dist);
                if (CheckPutVertex(cur_pair.vertex, cur_pair.edge, new_dist)) {
                    TRACE("CheckPutVertex returned true and new entry is added");
                    queue.push(element_t<Graph, distance_t>(new_dist, cur_pair.vertex,
                                    cur_vertex, cur_pair.edge));
                }
            }
            TRACE("Checking new neighbour of vertex " << graph_.str(cur_vertex) << " finished");
        }
        TRACE("All neighbours of vertex " << graph_.str(cur_vertex) << " processed");
    }

public:
    Dijkstra(const Graph &graph, DijkstraSettings settings, size_t max_vertex_number = size_t(-1)) :
        graph_(graph),
        settings_(settings),
        max_vertex_number_(max_vertex_number),
        finished_(false),
        vertex_number_(0),
        vertex_limit_exceeded_(false) {}

    Dijkstra(Dijkstra&& /*other*/) = default; 

    Dijkstra& operator=(Dijkstra&& /*other*/) = default;

    Dijkstra(const Dijkstra& /*other*/) = delete; 

    Dijkstra& operator=(const Dijkstra& /*other*/) = delete;

    bool finished() const {
        return finished_;
    }

    bool DistanceCounted(VertexId vertex) const {
        return distances_.find(vertex) != distances_.end();
    }

    distance_t GetDistance(VertexId vertex) const {
        VERIFY(DistanceCounted(vertex));
        return distances_.find(vertex)->second;
    }

    std::pair<distances_map_ci, distances_map_ci> GetDistances() const {
        distances_map_ci begin = distances_.begin();
        distances_map_ci end = distances_.end();
        return make_pair(begin, end);
    }

    void Run(VertexId start) {
        TRACE("Starting dijkstra run from vertex " << graph_.str(start));
        queue_t queue;
        Init(start, queue);
        TRACE("Priority queue initialized. Starting search");

        while (!queue.empty() && !finished()) {
            TRACE("Dijkstra iteration started");
            const element_t<Graph, distance_t>& next = queue.top();
            distance_t distance = next.distance;
            VertexId vertex = next.curr_vertex;

            prev_vert_map_[vertex] = std::pair<VertexId, EdgeId>(next.prev_vertex, next.edge_between);
            queue.pop();
            TRACE("Vertex " << graph_.str(vertex) << " with distance " << distance << " fetched from queue");

            if (DistanceCounted(vertex)) {
                TRACE("Distance to vertex " << graph_.str(vertex) << " already counted. Proceeding to next queue entry.");
                continue;
            }
            distances_.insert(make_pair(vertex, distance));

            TRACE("Vertex " << graph_.str(vertex) << " is found to be at distance "
                    << distance << " from vertex " << graph_.str(start));
            if (!CheckProcessVertex(vertex, distance)) {
                TRACE("Check for processing vertex failed. Proceeding to the next queue entry.");
                continue;
            }
            processed_vertices_.insert(vertex);
            AddNeighboursToQueue(vertex, distance, queue);
        }
        set_finished(true);
        TRACE("Finished dijkstra run from vertex " << graph_.str(start));
    }

    std::vector<EdgeId> GetShortestPathTo(VertexId vertex) {
        std::vector<EdgeId> path;
        if (prev_vert_map_.find(vertex) == prev_vert_map_.end())
            return path;

        VertexId curr_vertex = vertex;
        VertexId prev_vertex = get(prev_vert_map_, vertex).first;
        EdgeId edge = get(prev_vert_map_, curr_vertex).second;

        while (prev_vertex != VertexId(0)) {
            if (graph_.EdgeStart(edge) == prev_vertex)
                path.insert(path.begin(), edge);
            else
                path.push_back(edge);
            curr_vertex = prev_vertex;
            const auto& prev_v_e = get(prev_vert_map_, curr_vertex);
            prev_vertex = prev_v_e.first;
            edge = prev_v_e.second;
        }
        return path;
    }

    vector<VertexId> ReachedVertices() const {
        vector<VertexId> result;
        for (auto it = distances_.begin(); it != distances_.end(); ++it) {
            result.push_back(it->first);
        }
        return result;
    }

    const set<VertexId>& ProcessedVertices() const {
        return processed_vertices_;
    }

    bool VertexLimitExceeded() const {
        return vertex_limit_exceeded_;
    }

private:
    DECL_LOGGER("Dijkstra");
};

template<class Graph>
class DistanceCounter {
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef ComposedDijkstraSettings<Graph,
          LengthCalculator<Graph>,
              VertexProcessChecker<Graph>,
              VertexPutChecker<Graph>,
              ForwardNeighbourIteratorFactory<Graph>>  BaseDijkstraSettings;

public:
  DistanceCounter(const Graph& graph) :
    graph_(graph),
    dijkstra_(graph, BaseDijkstraSettings(
            LengthCalculator<Graph>(),
            VertexProcessChecker<Graph>(),
            VertexPutChecker<Graph>(),
            ForwardNeighbourIteratorFactory<Graph>())),
    ready_(false) {
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
  Dijkstra<Graph, BaseDijkstraSettings> dijkstra_;
  VertexId prev_;
  bool ready_;
};

}

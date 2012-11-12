//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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

	bool operator()(std::pair<distance_t, T> a, std::pair<distance_t, T> b) {
		if(a.first != b.first)
			return b.first < a.first;
		else
			return b.second < a.second;
	}
};

template<class Graph, typename distance_t = size_t>
class Dijkstra {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef std::map<VertexId, distance_t> distances_map;
	//	typedef map<VertexId, distance_t>::iterator distances_map_iterator;
	typedef typename distances_map::const_iterator distances_map_ci;

	const Graph &graph_;
	bool finished_;
	distances_map distances_;
	std::set<VertexId> processed_vertices_;

protected:
	const Graph& graph() {
		return graph_;
	}

public:

	bool finished() {
		return finished_;
	}

	Dijkstra(const Graph &graph) :
			graph_(graph), finished_(false) {
	}

	virtual ~Dijkstra() {
	}

	bool DistanceCounted(VertexId vertex) {
		return distances_.find(vertex) != distances_.end();
	}

	virtual void init(VertexId start) {
	}

	distance_t GetDistance(VertexId vertex) {
		VERIFY(DistanceCounted(vertex));
		return distances_.find(vertex)->second;
	}

	std::pair<distances_map_ci, distances_map_ci> GetDistances() {
		distances_map_ci begin = distances_.begin();
		distances_map_ci end = distances_.end();
		return make_pair(begin, end);
	}

	void set_finished(bool state) {
		finished_ = state;
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) {
		return true;
	}

	virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
		return true;
	}

	virtual distance_t GetLength(EdgeId edge) {
		return graph_.length(edge);
	}

	virtual std::vector<std::pair<VertexId, EdgeId> > Neighbours(VertexId vertex) {
		std::vector<std::pair<VertexId, EdgeId> > result;
    for (auto I = graph_.out_begin(vertex), E = graph_.out_end(vertex); I != E; ++I) {
      EdgeId edge = *I;
      result.push_back(make_pair(graph_.EdgeEnd(edge), edge));
		}
		return result;
	}

	void run(VertexId start) {
		TRACE("Starting dijkstra run from vertex " << graph_.str(start));
		TRACE("Initializing dijkstra priority queue");
		set_finished(false);
		distances_.clear();
		processed_vertices_.clear();
		init(start);
		std::priority_queue<std::pair<distance_t, VertexId> , std::vector<std::pair<distance_t,
				VertexId> > , ReverseDistanceComparator<distance_t, VertexId> > q;
		q.push(make_pair(0, start));
		TRACE("Priority queue initialized. Starting search");

		while (!q.empty() && !finished()) {
			TRACE("Dijkstra iteration started");
			auto next = q.top();
			q.pop();
			distance_t distance = next.first;
			VertexId vertex = next.second;
			TRACE(
					"Vertex " << graph_.str(vertex) << " with distance " << distance
							<< " fetched from queue");
			if (DistanceCounted(vertex)) {
				TRACE(
						"Distance to vertex " << graph_.str(vertex)
								<< " already counted. Proceeding to next queue entry.");
				continue;
			}
			distances_.insert(make_pair(vertex, distance));
			TRACE(
					"Vertex " << graph_.str(vertex) << " is found to be at distance "
							<< distance << " from vertex " << graph_.str(start));


			if (!CheckProcessVertex(vertex, distance)) {
				TRACE("Check for processing vertex failed. Proceeding to the next queue entry.");
				continue;
			}

			processed_vertices_.insert(vertex);
			auto neighbours = Neighbours(vertex);
			TRACE(
					"Neighbours of vertex " << graph_.str(vertex)
							<< " found. Iterating through neighbours and adding them to queue.");
			for (size_t i = 0; i < neighbours.size(); i++) {
				TRACE("Checking " << i << "th neighbour of vertex " << graph_.str(vertex) << " started");
				auto neighbour = neighbours[i];
//				TRACE("Which is " << neighbours[i]);
				if (!DistanceCounted(neighbour.first)) {
					TRACE("Adding new entry to queue");
					distance_t new_distance = GetLength(neighbour.second)
							+ distance;
					TRACE("Entry: vertex " << graph_.str(vertex) << " distance " << new_distance);
					if (CheckPutVertex(neighbour.first, neighbour.second,
							new_distance)) {
						TRACE("CheckPutVertex returned true and new entry is added");
						q.push(make_pair(new_distance, neighbour.first));
					}
				}
				TRACE("Checking " << i << "th neighbour of vertex " << graph_.str(vertex) << " finished");
			}
			TRACE("All neighbours of vertex " << graph_.str(vertex) << " processed");
		}
		set_finished(true);
		TRACE("Finished dijkstra run from vertex " << graph_.str(start));
	}

	std::vector<VertexId> ReachedVertices() {
		std::vector<VertexId> result;
		for (auto it = distances_.begin(); it != distances_.end(); ++it) {
			result.push_back(it->first);
		}
		return result;
	}

	const std::set<VertexId>& ProcessedVertices() {
		return processed_vertices_;
	}

private:
	DECL_LOGGER("Dijkstra");
};

template<class Graph>
class DistanceCounter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const Graph& graph_;
	Dijkstra<Graph> dijkstra_;
	VertexId prev_;
	bool ready_;
private:
	void EnsureFrom(VertexId from) {
		if (!ready_ || prev_ != from) {
			dijkstra_.run(from);
			ready_ = true;
			prev_ = from;
		}
	}

public:
	DistanceCounter(Graph &graph) :
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
};

template<class Graph, typename distance_t = size_t>
class UnorientedDijkstra: public Dijkstra<Graph, distance_t> {
private:
	typedef Dijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	UnorientedDijkstra(const Graph &graph) :
		super(graph) {
	}

	virtual ~UnorientedDijkstra() {
	}

	virtual std::vector<std::pair<VertexId, EdgeId> > Neighbours(VertexId vertex) {
		std::vector <std::pair<VertexId, EdgeId> > result;
		const Graph &g = this->graph();

    for (auto I = g.out_begin(vertex), E = g.out_end(vertex); I != E; ++I) {
      EdgeId edge = *I;
			result.push_back(make_pair(g.EdgeEnd(edge), edge));
    }

		std::vector <EdgeId> edges = g.IncomingEdges(vertex);
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(g.EdgeStart(edges[i]), edges[i]));
		}

    return result;
	}
};

template<class Graph, typename distance_t = size_t>
class BackwardDijkstra: public Dijkstra<Graph, distance_t> {
private:
	typedef Dijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	BackwardDijkstra(const Graph &graph) :
		super(graph) {
	}

	virtual ~BackwardDijkstra() {
	}

	virtual std::vector<std::pair<VertexId, EdgeId> > Neighbours(VertexId vertex) {
		TRACE("Starting to collect incoming edges for vertex " << this->graph().str(vertex));
		std::vector<std::pair<VertexId, EdgeId> > result;
		const Graph &g = this->graph();
		std::vector<EdgeId> edges = g.IncomingEdges(vertex);
		TRACE("Vector of incoming edges fetched from graph");
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(g.EdgeStart(edges[i]), edges[i]));
		}
		TRACE("Incoming edges info for vertex " << this->graph().str(vertex) << " constructed");
		return result;
	}
private:
	DECL_LOGGER("BackwardDijkstra");
};

template<class Graph, typename distance_t = size_t>
class BoundedDijkstra: public Dijkstra<Graph, distance_t> {
private:
	typedef Dijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	distance_t bound_;

public:
	BoundedDijkstra(const Graph &graph, distance_t bound) :
		super(graph), bound_(bound) {
	}

	virtual ~BoundedDijkstra() {
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) {
		if (length > bound_)
			return false;
		return true;
	}

	virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
		if (distance > bound_)
			return false;
		return true;
	}

};

}
#endif /* DIJKSTRA_HPP_ */

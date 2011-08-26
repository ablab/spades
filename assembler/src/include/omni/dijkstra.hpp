#ifndef DIJKSTRA_HPP_
#define DIJKSTRA_HPP_

#include <queue>

namespace omnigraph {

template<typename value_type>
class ReverseComparator {
public:
	bool operator()(value_type a, value_type b) {
		return a > b;
	}
};

template<class Graph, typename distance_t = size_t>
class Dijkstra {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<VertexId, distance_t> distances_map;
	//	typedef map<VertexId, distance_t>::iterator distances_map_iterator;
	typedef typename distances_map::const_iterator distances_map_ci;

	const Graph &graph_;
	bool finished_;
	map<VertexId, distance_t> distances_;
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
		assert(DistanceCounted(vertex));
		return distances_[vertex];
	}

	pair<distances_map_ci, distances_map_ci> GetDistances() {
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

	virtual vector<pair<VertexId, EdgeId>> Neighbours(VertexId vertex) {
		vector<pair<VertexId, EdgeId>> result;
		vector<EdgeId> edges = graph_.OutgoingEdges(vertex);
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(graph_.EdgeEnd(edges[i]), edges[i]));
		}
		return result;
	}

	void run(VertexId start) {
		set_finished(false);
		distances_.clear();
		init(start);
		priority_queue<pair<distance_t, VertexId> , vector<pair<distance_t,
				VertexId>> , ReverseComparator<pair<distance_t, VertexId>> > q;
		q.push(make_pair(0, start));
		while (!q.empty() && !finished()) {
			auto next = q.top();
			q.pop();
			distance_t distance = next.first;
			VertexId vertex = next.second;
			if (DistanceCounted(vertex)) {
				continue;
			}
			distances_.insert(make_pair(vertex, distance));
			if (!CheckProcessVertex(vertex, distance)) {
				continue;
			}
			auto neighbours = Neighbours(vertex);
			for (size_t i = 0; i < neighbours.size(); i++) {
				auto neighbour = neighbours[i];
				if (!DistanceCounted(neighbour.first)) {
					distance_t new_distance = GetLength(neighbour.second)
							+ distance;
					if (CheckPutVertex(neighbour.first, neighbour.second,
							new_distance)) {
						q.push(make_pair(new_distance, neighbour.first));
					}
				}
			}
		}
		set_finished(true);
	}

	vector<VertexId> VisitedVertices() {
		vector<VertexId> result;
		for (auto it = distances_.begin(); it != distances_.end(); ++it) {
			result.push_back(it->first);
		}
		return result;
	}

};

template<class Graph>
class DistanceCounter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Graph& graph_;
	Dijkstra<Graph> dijkstra_;
	VertexId prev_;
	bool ready_;
private:
	void EnsureFrom(VertexId from) {
		if(!ready_ || prev_ != from) {
			dijkstra_.run(from);
			ready_ = true;
			prev_ = from;
		}
	}

public:
	DistanceCounter(Graph &graph) :graph_(graph), dijkstra_(graph), ready_(false){
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

	virtual vector<pair<VertexId, EdgeId>> Neighbours(VertexId vertex) {
		vector < pair < VertexId, EdgeId >> result;
		const Graph &g = this->graph();
		vector < EdgeId > edges = g.OutgoingEdges(vertex);
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(g.EdgeEnd(edges[i]), edges[i]));
		}
		edges = g.IncomingEdges(vertex);
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

	virtual vector<pair<VertexId, EdgeId>> Neighbours(VertexId vertex) {
		vector <pair<VertexId, EdgeId>> result;
		const Graph &g = this->graph();
		vector<EdgeId> edges = g.OutgoingEdges(vertex);
		edges = g.IncomingEdges(vertex);
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(g.EdgeStart(edges[i]), edges[i]));
		}
		return result;
	}
};

template<class Graph, typename distance_t = size_t>
class BoundedDijkstra: public Dijkstra<Graph, distance_t> {
private:
	typedef Dijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	distance_t bound_;

public:
	BoundedDijkstra(Graph &graph, distance_t bound) :
		super(graph), bound_(bound) {
	}

	virtual ~BoundedDijkstra() {
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) {
		if (length > bound_) return false;
		return true;
	}

	virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
		if (distance > bound_) return false;
		return true;
	}

};
}
#endif /* DIJKSTRA_HPP_ */

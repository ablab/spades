#ifndef DIJKSTRA_HPP_
#define DIJKSTRA_HPP_

#include <queue>

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

	Graph &graph_;
	bool finished_;
	map<VertexId, distance_t> distances_;
public:

	Graph &GetGraph() {
		return graph_;
	}

	bool finished() {
		return finished_;
	}

	Dijkstra(Graph &graph) :
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

template<class Graph, typename distance_t = size_t>
class UnorientedDijkstra: public Dijkstra<Graph, distance_t> {
private:
	typedef Dijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	UnorientedDijkstra(Graph &graph) :
		super(graph) {
	}

	virtual ~UnorientedDijkstra() {
	}

	virtual vector<pair<VertexId, EdgeId>> Neighbours(VertexId vertex) {
		vector < pair < VertexId, EdgeId >> result;
		Graph &g = this->GetGraph();
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
#endif /* DIJKSTRA_HPP_ */

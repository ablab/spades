#ifndef DIJKSTRA_HPP_
#define DIJKSTRA_HPP_

#include <queue>
#include <map>
#include <set>

namespace omnigraph {

template<typename distance_t, typename T, class Comparator>
class ReverseDistanceComparator {
private:
	Comparator comparator_;
public:
	ReverseDistanceComparator(Comparator comparator): comparator_(comparator) {
	}

	bool operator()(std::pair<distance_t, T> a, std::pair<distance_t, T> b) {
		if(a.first != b.first)
			return b.first < a.first;
		else
			return comparator_(b.second, a.second);
	}
};

template<class Graph, typename distance_t = size_t>
class Dijkstra {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef std::map<VertexId, distance_t, typename Graph::Comparator> distances_map;
	//	typedef map<VertexId, distance_t>::iterator distances_map_iterator;
	typedef typename distances_map::const_iterator distances_map_ci;

	const Graph &graph_;
	bool finished_;
	distances_map distances_;
	std::set<VertexId, typename Graph::Comparator> processed_vertices_;

protected:
	const Graph& graph() {
		return graph_;
	}

public:

	bool finished() {
		return finished_;
	}

	Dijkstra(const Graph &graph) :
			graph_(graph), finished_(false), distances_(
					graph_.ReliableComparatorInstance()), processed_vertices_(
					graph_.ReliableComparatorInstance()) {
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
		return distances_[vertex];
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
		std::vector<EdgeId> edges = graph_.OutgoingEdges(vertex);
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(graph_.EdgeEnd(edges[i]), edges[i]));
		}
		return result;
	}

	void run(VertexId start) {
		TRACE("Starting dijkstra run from vertex " << start);
		TRACE("Initializing dijkstra priority queue");
		set_finished(false);
		distances_.clear();
		processed_vertices_.clear();
		init(start);
		std::priority_queue<std::pair<distance_t, VertexId> , std::vector<std::pair<distance_t,
				VertexId> > , ReverseDistanceComparator<distance_t, VertexId, typename Graph::Comparator> > q(graph_.ReliableComparatorInstance());
		q.push(make_pair(0, start));
		TRACE("Priority queue initialized. Starting search");

		while (!q.empty() && !finished()) {
			TRACE("Dijkstra iteration started");
			auto next = q.top();
			q.pop();
			distance_t distance = next.first;
			VertexId vertex = next.second;
			TRACE(
					"Vertex " << vertex << " with distance " << distance
							<< "fetched from queue");
			if (DistanceCounted(vertex)) {
				TRACE(
						"Distance to vertex " << vertex
								<< " already counted. Proceeding to next queue entry.");
				continue;
			}
			distances_.insert(make_pair(vertex, distance));
			TRACE(
					"Vertex " << vertex << " is found to be at distance "
							<< distance << " from vertex " << start);


			if (!CheckProcessVertex(vertex, distance)) {
				TRACE("Check for processing vertex failed. Proceeding to the next queue entry.");
				continue;
			} else {
				processed_vertices_.insert(vertex);
			}

			auto neighbours = Neighbours(vertex);
			TRACE(
					"Neighbours of vertex " << vertex
							<< " found. Iterating through neighbours and adding them to queue.");
			for (size_t i = 0; i < neighbours.size(); i++) {
				TRACE("Checking " << i << "th neighbour of vertex " << vertex << " started");
				auto neighbour = neighbours[i];
				TRACE("Which is " << neighbours[i]);
				if (!DistanceCounted(neighbour.first)) {
					TRACE("Adding new entry to queue");
					distance_t new_distance = GetLength(neighbour.second)
							+ distance;
					TRACE("Entry: vertex " << vertex << " distance " << new_distance);
					if (CheckPutVertex(neighbour.first, neighbour.second,
							new_distance)) {
						TRACE("CheckPutVertex returned true and new entry is added");
						q.push(make_pair(new_distance, neighbour.first));
					}
				}
				TRACE("Checking " << i << "th neighbour of vertex " << vertex << " finished");
			}
			TRACE("All neighbours of vertex " << vertex << " processed");
		}
		set_finished(true);
		TRACE("Finished dijkstra run from vertex " << start);
	}

	std::vector<VertexId> ReachedVertices() {
		std::vector<VertexId> result;
		for (auto it = distances_.begin(); it != distances_.end(); ++it) {
			result.push_back(it->first);
		}
		return result;
	}

	const std::set<VertexId, typename Graph::Comparator>& ProcessedVertices() {
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
		std::vector <EdgeId> edges = g.OutgoingEdges(vertex);
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

	virtual std::vector<std::pair<VertexId, EdgeId> > Neighbours(VertexId vertex) {
		TRACE("Starting to collect incoming edges for vertex " << vertex);
		std::vector<std::pair<VertexId, EdgeId> > result;
		const Graph &g = this->graph();
		std::vector<EdgeId> edges = g.IncomingEdges(vertex);
		TRACE("Vector of incoming edges fetched from graph");
		for (size_t i = 0; i < edges.size(); i++) {
			result.push_back(make_pair(g.EdgeStart(edges[i]), edges[i]));
		}
		TRACE("Incoming edges info for vertex " << vertex << " constructed");
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

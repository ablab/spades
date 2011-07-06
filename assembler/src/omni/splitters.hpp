#ifndef SPLITTERS_HPP_
#define SPLITTERS_HPP_

#include "dijkstra.hpp"
namespace omnigraph {
template<class Element>
class GraphSplitter {
public:
	virtual vector<Element> NextComponent() = 0;

	virtual bool Finished() = 0;

	virtual ~GraphSplitter() {

	}
};

template<class Graph>
class ComponentFinder: public UnorientedDijkstra<Graph, size_t> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef UnorientedDijkstra<Graph, size_t> super;
	set<EdgeId> &edges_;

public:
	ComponentFinder(Graph &g, set<EdgeId> &edges) :
		super(g), edges_(edges) {
	}

	virtual ~ComponentFinder() {
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, size_t length) {
		return edges_.count(edge) != 0;
	}
};

template<class Graph>
class NeighbourhoodFinder: public UnorientedDijkstra<Graph, size_t> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef UnorientedDijkstra<Graph, size_t> super;
	set<EdgeId> &edges_;
	const size_t bound_;

public:
	NeighbourhoodFinder(Graph &g, set<EdgeId> &edges, size_t bound) :
		super(g), edges_(edges), bound_(bound) {
	}

	virtual ~NeighbourhoodFinder() {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		return distance <= bound_;
	}

	virtual size_t GetLength(EdgeId edge) {
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
	typedef UnorientedDijkstra<Graph, size_t> super;
	const set<VertexId> &subgraph_;

public:
	SubgraphDijkstra(Graph &g, const set<VertexId> &subgraph) :
		super(g), subgraph_(subgraph) {
	}

	virtual ~SubgraphDijkstra() {
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, size_t length) {
		return subgraph_.count(vertex) != 0;
	}

};

template<class Graph>
class ErrorComponentSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
	set<EdgeId> black_edges_;
	SmartEdgeIterator<omnigraph::ObservableGraph<VertexId, EdgeId> > iterator_;
	set<VertexId> visited_;

public:
	ErrorComponentSplitter(Graph &graph, const set<EdgeId> &black_edges) :
		graph_(graph), black_edges_(black_edges),
				iterator_(graph.SmartEdgeBegin()) {
	}

	virtual ~ErrorComponentSplitter() {
	}

	set<VertexId> FindComponent(VertexId start_vertex) {
		ComponentFinder<Graph> cf(graph_, black_edges_);
		cf.run(start_vertex);
		vector < VertexId > result = cf.VisitedVertices();
		return set<VertexId> (result.begin(), result.end());
	}

	set<VertexId> FindNeighbourhood(VertexId start, size_t bound) {
		NeighbourhoodFinder<Graph> nf(graph_, black_edges_, bound);
		nf.run(start);
		vector < VertexId > result = nf.VisitedVertices();
		return set<VertexId> (result.begin(), result.end());
	}

	size_t FindDiameter(const set<VertexId> &component) {
		size_t result = 0;
		SubgraphDijkstra<Graph> sd(graph_, component);
		for (auto it = component.begin(); it != component.end(); ++it) {
			sd.run(*it);
			auto bounds = sd.GetDistances();
			for (auto it = bounds.first; it != bounds.second; ++it) {
				result = std::max(result, it->second);
			}
		}
		return result;
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			assert(false);
			return vector<VertexId> ();
		}
		EdgeId next = *iterator_;
		++iterator_;
		set < VertexId > component = FindComponent(graph_.EdgeEnd(next));
		size_t component_size = FindDiameter(component);
		set < VertexId > neighbourhood = FindNeighbourhood(
				graph_.EdgeEnd(next), 1.5 * component_size);
		visited_.insert(component.begin(), component.end());
		return vector<VertexId> (neighbourhood.begin(), neighbourhood.end());
	}

	virtual bool Finished() {
		while (!iterator_.IsEnd()) {
			if (black_edges_.find(*iterator_) != black_edges_.end()
					&& visited_.find(graph_.EdgeEnd(*iterator_))
							== visited_.end()) {
				return false;
			}
			++iterator_;
		}
		return true;
	}

};

template<class Graph>
class ShortEdgeComponentFinder: public UnorientedDijkstra<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	size_t bound_;
public:
	ShortEdgeComponentFinder(Graph &graph, size_t bound) :
		UnorientedDijkstra<Graph> (graph), bound_(bound) {
	}

	virtual ~ShortEdgeComponentFinder() {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		return distance == 0;
	}

	virtual size_t GetLength(EdgeId edge) {
		if (this->graph().length(edge) <= bound_)
			return 0;
		else
			return 1;
	}
};

template<class Graph>
class LongEdgesSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
	SmartVertexIterator<omnigraph::ObservableGraph<VertexId, EdgeId> >
			iterator_;
	set<VertexId> visited_;
	size_t bound_;

public:
	LongEdgesSplitter(Graph &graph, size_t bound) :
		graph_(graph), bound_(bound), iterator_(graph.SmartVertexBegin()) {
	}

	virtual ~LongEdgesSplitter() {
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			assert(false);
			return vector<VertexId> ();
		}
		VertexId next = *iterator_;
		ShortEdgeComponentFinder<Graph> cf(graph_, bound_);
		cf.run(next);
		vector<VertexId> result;
		for(auto it = cf.begin(); it != cf.end(); ++it) {
			result.push_back(it->first);
			if(it->second == 0) {
				iterator_.erase(it->first);
			}
		}
		return result;
	}

	virtual bool Finished() {
		return iterator_.IsEnd();
	}

};

}

#endif /* SPLITTERS_HPP_ */

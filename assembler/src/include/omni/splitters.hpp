#ifndef SPLITTERS_HPP_
#define SPLITTERS_HPP_

#include "dijkstra.hpp"
namespace omnigraph {
template<class Element>
class GraphSplitter {
public:
	virtual vector<Element> NextComponent() = 0;

	virtual bool Finished() = 0;

	virtual string ComponentName() {
		return "";
	}

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
	ComponentFinder(const Graph &g, set<EdgeId> &edges) :
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
	NeighbourhoodFinder(const Graph &g, set<EdgeId> &edges, size_t bound) :
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
	SubgraphDijkstra(const Graph &g, const set<VertexId> &subgraph) :
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

	const Graph &graph_;
	set<EdgeId> black_edges_;
	typename Graph::SmartEdgeIt iterator_;
	set<VertexId> visited_;

public:
	ErrorComponentSplitter(const Graph &graph, const set<EdgeId> &black_edges) :
		graph_(graph), black_edges_(black_edges),
				iterator_(graph.SmartEdgeBegin()) {
		TRACE("ErrorComponentSplitter created and SmartIterator initialized");
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
		VertexId current = *(component.begin());
		for (size_t i = 0; i < 4; i++) {
			pair < VertexId, size_t > next = GetFarthest(current, component);
			current = next.first;
			result = next.second;
		}
		return result;
	}

	pair<VertexId, size_t> GetFarthest(VertexId v,
			const set<VertexId> &component) {
		SubgraphDijkstra<Graph> sd(graph_, component);
		sd.run(v);
		pair < VertexId, size_t > result(v, 0);
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
			return vector<VertexId> ();
		}
		EdgeId next = *iterator_;
		++iterator_;
		set < VertexId > component = FindComponent(graph_.EdgeEnd(next));
		TRACE(
				"Error edges component constructed. It contains "
						<< component.size() << " vertices");
		size_t component_size = FindDiameter(component);
		TRACE("Diameter of component is " << component_size);
		set < VertexId > neighbourhood = FindNeighbourhood(
				graph_.EdgeEnd(next), 1.5 * component_size);
		TRACE(
				"Error edges component neighborhood constructed. It contains "
						<< neighbourhood.size() << " vertices");
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
class ShortEdgeComponentNeighbourhoodFinder: public UnorientedDijkstra<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	size_t bound_;
public:
	ShortEdgeComponentNeighbourhoodFinder(const Graph &graph, size_t bound) :
		UnorientedDijkstra<Graph> (graph), bound_(bound) {
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
class LongEdgesInclusiveSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph &graph_;
	erasable_priority_queue<VertexId> queue_;
	//	SmartVertexIterator<omnigraph::ObservableGraph<VertexId, EdgeId> >
	//			iterator_;
	set<VertexId> visited_;
	size_t bound_;

public:
	LongEdgesInclusiveSplitter(const Graph &graph, size_t bound) :
		graph_(graph), queue_(graph.begin(), graph.end()), /*iterator_(graph.SmartVertexBegin()), */
		bound_(bound) {
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			VERIFY(false);
			return vector<VertexId> ();
		}
		VertexId next = queue_.top();
		TRACE("Search started");
		queue_.pop();
		ShortEdgeComponentNeighbourhoodFinder<Graph> cf(graph_, bound_);
		cf.run(next);
		TRACE("Search finished");
		vector < VertexId > result = cf.VisitedVertices();
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
	typedef UnorientedDijkstra<Graph, distance_t> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	static const size_t inf = 100000;

	size_t max_size_;

	size_t edge_length_bound_;

	size_t current_;

public:
	CountingDijkstra(const Graph &graph, size_t max_size,
			size_t edge_length_bound) :
		super(graph), max_size_(max_size),
				edge_length_bound_(edge_length_bound), current_(0) {
	}

	virtual ~CountingDijkstra() {
	}

	virtual bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) {
		if (current_ < max_size_) {
			current_++;
		}
		if (current_ < max_size_ || length >= inf) {
			return true;
		}
		return false;
	}

	virtual bool CheckProcessVertex(VertexId vertex, distance_t distance) {
		return distance < inf && current_ < max_size_;
	}

	virtual size_t GetLength(EdgeId edge) {
		if (this->graph().length(edge) <= edge_length_bound_)
			return this->graph().length(edge);
		else
			return inf;
	}
};

template<class Graph>
class ShortEdgeComponentFinder: public UnorientedDijkstra<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	size_t bound_;
public:
	ShortEdgeComponentFinder(const Graph &graph, size_t bound) :
		UnorientedDijkstra<Graph> (graph), bound_(bound) {
	}

	virtual bool CheckPutVertex(VertexId vertex, size_t distance) {
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
class ReliableSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph &graph_;
	set<VertexId> visited_;
	size_t max_size_;
	size_t edge_length_bound_;
	typename Graph::VertexIterator current_;

	void SkipVisited() {
		while (current_ != graph_.end() && visited_.count(*current_) == 1) {
			++current_;
		}
	}

public:
	ReliableSplitter(const Graph &graph, size_t max_size,
			size_t edge_length_bound) :
		graph_(graph), max_size_(max_size),
				edge_length_bound_(edge_length_bound), current_(graph.begin()) {
		TRACE(
				"Long edges splitter created and queue filled with all graph vertices");
	}

	virtual ~ReliableSplitter() {
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			VERIFY(false);
			return vector<VertexId> ();
		}
		TRACE("Search started");
		CountingDijkstra<Graph> cf(graph_, max_size_, edge_length_bound_);
		cf.run(*current_);
		TRACE("Search finished");
		vector < VertexId > result = cf.VisitedVertices();
		visited_.insert(result.begin(), result.end());
		TRACE("Component vector filled");
		SkipVisited();
		return result;
	}

	virtual bool Finished() {
		return current_ == graph_.end();
	}
};

template<class Graph>
class ReliableSplitterAlongGenome: public GraphSplitter<
		typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph &graph_;

	size_t max_size_;
	size_t edge_length_bound_;
	set<VertexId> last_component_;
	size_t current_index_;
	MappingPath<EdgeId> genome_path_;
	Range covered_range_;
	bool start_processed_;

	bool EdgeCovered(EdgeId edge) {
		return last_component_.count(
				graph_.EdgeStart(genome_path_[current_index_].first)) == 1
				&& last_component_.count(
						graph_.EdgeEnd(genome_path_[current_index_].first)) == 1;
	}

	void SkipVisited() {
		covered_range_.start_pos
				= genome_path_[current_index_].second.initial_range.start_pos;
		covered_range_.end_pos
				= genome_path_[current_index_].second.initial_range.end_pos;
		while (current_index_ != genome_path_.size() && EdgeCovered(
				genome_path_[current_index_].first)) {
			covered_range_.end_pos
					= genome_path_[current_index_].second.initial_range.end_pos;
			++current_index_;
		}
	}

public:
	string ComponentName() {
		stringstream ss;
		ss << covered_range_.start_pos << "_" << covered_range_.end_pos;
		return ss.str();
	}

	ReliableSplitterAlongGenome(const Graph &graph, size_t max_size,
			size_t edge_length_bound, MappingPath<EdgeId> genome_path) :
		graph_(graph), max_size_(max_size),
				edge_length_bound_(edge_length_bound), current_index_(0),
				genome_path_(genome_path), covered_range_(0, 0), start_processed_(false) {

	}

	virtual ~ReliableSplitterAlongGenome() {
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			VERIFY(false);
			return vector<VertexId> ();
		}
		TRACE("Search started");
		CountingDijkstra<Graph> cf(graph_, max_size_, edge_length_bound_);
		if(start_processed_)
			cf.run(graph_.EdgeEnd(genome_path_[current_index_].first));
		else {
			cf.run(graph_.EdgeStart(genome_path_[current_index_].first));
			start_processed_ = true;
		}
		TRACE("Search finished");
		vector < VertexId > result = cf.VisitedVertices();
		last_component_.clear();
		last_component_.insert(result.begin(), result.end());
		VERIFY(EdgeCovered(genome_path_[current_index_].first));
		last_component_.insert(result.begin(), result.end());
		TRACE("Component vector filled");
		size_t prev_index = current_index_;
		SkipVisited();
		if(prev_index + 1 != current_index_) {
			start_processed_ = true;
		} else if(!start_processed_) {
			current_index_ = prev_index;
			start_processed_ = true;
		} else {
			start_processed_ = false;
		}
		return result;
	}

	virtual bool Finished() {
		return current_index_ == genome_path_.size();
	}
};

template<class Graph>
class LongEdgesExclusiveSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph &graph_;
	erasable_priority_queue<VertexId> queue_;
	//	SmartVertexIterator<omnigraph::ObservableGraph<VertexId, EdgeId> >
	//			iterator_;
	set<VertexId> visited_;
	size_t bound_;

public:
	LongEdgesExclusiveSplitter(const Graph &graph, size_t bound) :
		graph_(graph), queue_(graph.begin(), graph.end()), /*iterator_(graph.SmartVertexBegin()), */
		bound_(bound) {
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			VERIFY(false);
			return vector<VertexId> ();
		}
		VertexId next = queue_.top();
		queue_.pop();
		ShortEdgeComponentFinder<Graph> cf(graph_, bound_);
		cf.run(next);
		vector < VertexId > result = cf.VisitedVertices();
		for (auto it = result.begin(); it != result.end(); ++it) {
			queue_.erase(*it);
		}
		return result;
	}

	virtual bool Finished() {
		//		return iterator_.IsEnd();
		return queue_.empty();
	}

};

template<class Element>
class AbstractFilter {
public:
	virtual ~AbstractFilter() {
	}

	virtual bool Check(Element &element) const = 0;
};

template<class Graph>
class ComponentSizeFilter: public AbstractFilter<vector<
		typename Graph::VertexId>> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& graph_;
	size_t max_length_;

	//	bool CheckYellow() {
	//
	//	}
	//
public:
	ComponentSizeFilter(const Graph &graph, size_t max_length) :
		graph_(graph), max_length_(max_length) {
	}

	virtual ~ComponentSizeFilter() {
	}

	virtual bool Check(vector<VertexId> &vertices) const {
		if (vertices.size() <= 4)
			return false;
		set<VertexId> component(vertices.begin(), vertices.end());
		for (auto iterator = vertices.begin(); iterator != vertices.end(); ++iterator) {
			vector<EdgeId> edges = graph_.OutgoingEdges(*iterator);
			for (auto edge_iterator = edges.begin(); edge_iterator
					!= edges.end(); edge_iterator++) {
				if (component.count(graph_.EdgeEnd(*edge_iterator)) == 1
						&& graph_.length(*edge_iterator) <= max_length_) {
					return true;
				}
			}
		}
		return false;
	}
};

template<class Graph>
class FilteringSplitterWrapper: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::VertexId VertexId;
	GraphSplitter<typename Graph::VertexId> &inner_splitter_;
	vector<VertexId> next;
	const AbstractFilter<vector<VertexId>> &checker_;
	bool ready;
public:
	FilteringSplitterWrapper(
			GraphSplitter<typename Graph::VertexId> &inner_splitter,
			const AbstractFilter<vector<VertexId>> &checker) :
		inner_splitter_(inner_splitter), checker_(checker), ready(false) {
	}

	virtual ~FilteringSplitterWrapper() {
	}

	virtual string ComponentName() {
		return inner_splitter_.ComponentName();
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			VERIFY(false);
			return vector<VertexId> ();
		}
		ready = false;
		return next;
	}

	virtual bool Finished() {
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

};

}

#endif /* SPLITTERS_HPP_ */

#pragma once

namespace omnigraph {

using std::set;
using std::map;
using std::vector;
using std::pair;
using std::queue;

template<class Graph>
class FlowGraph {
public:
	typedef size_t VertexId;
	typedef pair<VertexId, VertexId> EdgeId;

private:
	typedef typename Graph::VertexId OuterVertexId;
	map<OuterVertexId, VertexId> vertex_mapping_;
	map<VertexId, map<VertexId, int>> capacities_;
	set<VertexId> vertices_;
	size_t vertex_number_;
	VertexId source_;
	VertexId sink_;

	VertexId AddVertex() {
		vertex_number_++;
		return vertex_number_ - 1;
	}

	void PushFlow(EdgeId edge, int capacity) {
		VERIFY(capacities_[EdgeStart(edge)][EdgeEnd(edge)] >= capacity);
		capacities_[EdgeStart(edge)][EdgeEnd(edge)] -= capacity;
		capacities_[EdgeEnd(edge)][EdgeStart(edge)] += capacity;
	}

	void AddEdge(VertexId first, VertexId second, int capacity = 10000) {
		capacities_[first][second] += capacity; // operator [] creates entry with default values in case argument is not in keyset
		capacities_[second][first] += 0;
	}

public:
	FlowGraph() :
			vertex_number_(0), source_(AddVertex()), sink_(AddVertex()) {
	}

	VertexId GetCorrespondingVertex(OuterVertexId v) const {
		return vertex_mapping_.find(v)->second;
	}

	bool HasCorrespondingVertex(OuterVertexId v) const {
		return vertex_mapping_.find(v) == vertex_mapping_.end();
	}

	VertexId AddVertex(OuterVertexId vertex) {
		VertexId new_vertex = AddVertex();
		vertex_mapping_[vertex] = new_vertex;
		return new_vertex;
	}

	void AddEdge(OuterVertexId outer_first, OuterVertexId outer_second,
			int capacity = 10000) {
		VERIFY(
				vertex_mapping_.find(outer_first) != vertex_mapping_.end()
						&& vertex_mapping_.find(outer_second)
								!= vertex_mapping_.end());
		VertexId first = vertex_mapping_[outer_first];
		VertexId second = vertex_mapping_[outer_second];
		AddEdge(first, second, capacity);
	}

	void AddSource(OuterVertexId vertex, int capacity) {
		AddEdge(source_, GetCorrespondingVertex(vertex), capacity = 10000);
	}

	void AddSink(OuterVertexId vertex, int capacity) {
		AddEdge(GetCorrespondingVertex(vertex), sink_, capacity);
	}

	VertexId Source() const {
		return source_;
	}

	VertexId Sink() const {
		return sink_;
	}

	bool Connected(VertexId start, VertexId end) const {
		return capacities_.find(start) != capacities_.end()
				&& capacities_.find(start)->second.find(end)
						!= capacities_.find(start)->second.end()
				&& capacities_.find(start)->second.find(end)->second > 0;
	}

	vector<EdgeId> OutgoingEdges(VertexId v) const {
		vector<EdgeId> result;
		const map<VertexId, int> &outgoing = capacities_.find(v)->second;
		for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
			if (it->second > 0) {
				result.push_back(make_pair(v, it->first));
			}
		}
		return result;
	}

	vector<EdgeId> IncomingEdges(VertexId v) const {
		vector<EdgeId> result;
		const map<VertexId, int> &outgoing = capacities_.find(v)->second;
		for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
			if (Connected(it->first, v)) {
				result.push_back(make_pair(it->first, v));
			}
		}
		return result;
	}

	size_t OutgoingEdgesCount(VertexId v) const {
		return OutgoingEdges(v).size();
	}

	size_t IncomingEdgesCount(VertexId v) const {
		return IncomingEdges(v).size();
	}

	VertexId EdgeStart(EdgeId edge) const {
		return edge.first;
	}

	VertexId EdgeEnd(EdgeId edge) const {
		return edge.second;
	}

	set<VertexId>::iterator begin() const {
		return vertices_.begin();
	}

	set<VertexId>::iterator end() const {
		return vertices_.end();
	}

	int GetCapacity(VertexId first, VertexId second) const {
		auto it1 = capacities_.find(first);
		if (it1 == capacities_.end())
			return 0;
		auto it2 = it1->second.find(second);
		if (it2 == it1->second.end())
			return 0;
		return it2->second;
	}

	void PushFlow(vector<VertexId> path, int capacity) {
		size_t n = path.size();
		VERIFY(path[0] == source_ && path[n - 1] == sink_);
		for (size_t i = 0; i + 1 < n; i++) {
			PushFlow(make_pair(path[i], path[i + 1]), capacity);
		}
	}
};

template<class Graph>
class BFS {
private:
	const Graph &graph_;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	vector<VertexId> RestoreAnswer(VertexId start, VertexId end, const map<VertexId, VertexId> &prev) {
		vector<VertexId> result;
		result.push_back(end);
		VertexId current = end;
		while(current != start) {
			current = prev.find(current)->second;
			result.push_back(current);
		}
		return vector<VertexId>(result.rbegin(), result.rend());
	}

public:
	BFS(const Graph &graph) :
			graph_(graph) {
	}

	vector<VertexId> Go(VertexId start, VertexId finish) {
		queue<VertexId> q;
		q.push(start);
		map<VertexId, VertexId> prev;
		prev[start] = start;
		while(!q.empty()) {
			VertexId current = q.front();
			q.pop();
			vector<EdgeId> outgoing = graph_.OutgoingEdges(current);
			for(auto it = outgoing.begin(); it != outgoing.end(); ++it) {
				if(prev.find(it->second) == prev.end()) {
					q.push(it->second);
					prev[it->second] = current;
				}
				if(it->second == finish) {
					return RestoreAnswer(start, finish, prev);
				}
			}
		}
		return vector<VertexId>();
	}
};

template<class Graph>
class MaxFlowFinder {
private:
	FlowGraph<Graph> &graph_;
	typedef typename FlowGraph<Graph>::VertexId VertexId;
	typedef typename FlowGraph<Graph>::EdgeId EdgeId;

	int MaxCapacity(vector<VertexId> path) {
		VERIFY(path.size() > 2)
		int result = graph_.GetCapacity(path[0], path[1]);
		for(size_t i = 1; i + 1 < path.size(); i++) {
			result = std::max(result, graph_.GetCapacity(i, i + 1));
		}
		return result;
	}

public:
	MaxFlowFinder(FlowGraph<Graph> &graph) :
			graph_(graph) {
	}

	void Find() {
		BFS<FlowGraph<Graph> > bfs(graph_);
		while(true) {
			vector<VertexId> path = bfs.Go(graph_.Source(), graph_.Sink());
			if(path.size() == 0)
				break;
			int capacity = MaxCapacity(path);
			VERIFY(capacity > 0)
			graph_.PushFlow(path, capacity);
		}
	}
};

template<class Graph>
class TopSorter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;

	void Find(VertexId v, vector<VertexId> &result, set<VertexId> &visited) {
		visited.insert(v);
		vector<EdgeId> outgoing = graph_.OutgoingEdges(v);
		for(auto it = outgoing.begin(); it != outgoing.end(); ++it) {
			VertexId next = graph_.EdgeEnd(*it);
			if(visited.count(next) == 0) {
				Find(next, result, visited);
			}
		}
		result.push_back(v);
	}

public:
	TopSorter(const Graph &graph) : graph_(graph) {
	}

	vector<VertexId> Sort() {
		vector<VertexId> result;
		set<VertexId> visited;
		for(auto it = graph_.begin(); it != graph_.end(); ++it) {
			if(visited.count(*it) == 0) {
				Find(*it, result, visited);
			}
		}
		return result;
	}
};


template<class Graph>
class ReverseDFSComponentFinder {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const Graph &graph_;

	void Find(VertexId v, map<VertexId, size_t> &result, size_t cc) {
		result[v] = cc;
		vector<EdgeId> incoming= graph_.IncomingEdges(v);
		for(auto it = incoming.begin(); it != incoming.end(); ++it) {
			VertexId next = graph_.EdgeStart(*it);
			if(result.count(next) == 0) {
				Find(next, result, cc);
			}
		}
	}
public:
	ReverseDFSComponentFinder(const Graph &graph) : graph_(graph) {
	}

	map<VertexId, size_t> Find(const vector<VertexId> &order) {
		size_t cc = 0;
		map<VertexId, size_t> result;
		for(auto it = order.rbegin(); it != order.rend(); ++it) {
			if(result.count(*it)) {
				Find(*it, result, cc);
				cc++;
			}
		}
		return result;
	}
};

template<class Graph>
class StroglyConnectedComponentFinder {
private:
	typedef typename Graph::VertexId VertexId;
	const Graph &graph_;
	bool ready_;
public:
	StroglyConnectedComponentFinder(const Graph &graph) :
			graph_(graph), ready_(false) {
	}

	map<VertexId, size_t> ColourComponents() {
		map<VertexId, size_t> result;
		vector<VertexId> order = TopSorter<Graph>(graph_).Sort();
		return ReverseDFSComponentFinder<Graph>(graph_).Find(order);
	}
};

template<class Graph>
class MaxFlowECRemover {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph& graph_;
	size_t max_length_;
	size_t uniqueness_length_;
	size_t plausibility_length_;
	EdgeRemover<Graph>& edge_remover_;

	bool IsTerminal(VertexId vertex) {
		return graph_.OutgoingEdgeCount(vertex)
				+ graph_.IncomingEdgeCount(vertex) == 1;
	}

	bool IsTip(EdgeId edge) {
		VertexId start = graph_.EdgeStart(edge);
		VertexId end = graph_.EdgeEnd(edge);
		return IsTerminal(start) || IsTerminal(end);
	}

	bool IsSuspicious(EdgeId edge) {
		return graph_.length(edge) <= max_length_ && !IsTip(edge);
	}

	set<EdgeId> CollectUnusedEdges(set<VertexId> component, FlowGraph<Graph> fg,
			const map<typename FlowGraph<Graph>::VertexId, size_t> &colouring) {
		set<EdgeId> result;
		for (auto it_start = component.begin(); it_start != component.end();
				++it_start) {
			VertexId start = *it_start;
			auto outgoing = graph_.OutgoingEdges(start);
			for (auto it_edge = outgoing.begin(); it_edge != outgoing.end();
					++it_edge) {
				EdgeId edge = *it_edge;
				VertexId end = graph_.EdgeEnd(edge);
				if (component.count(end) == 1 && IsSuspicious(edge)
						&& colouring.find(fg.GetCorrespondingVertex(start))->second
								!= colouring.find(
										fg.GetCorrespondingVertex(end))->second) {
					result.insert(edge);
				}
			}
		}
		return result;
	}

	bool CheckCompleteFlow(FlowGraph<Graph> &fg) {
		return fg.OutgoingEdges(fg.Source()).size() == 0
				&& fg.IncomingEdges(fg.Source()).size() == 0;
	}

	bool IsPlausible(EdgeId edge) {
		return graph_.length(edge) >= plausibility_length_ && !IsTip(edge);
	}

	bool IsUnique(EdgeId edge) {
		return graph_.length(edge) >= uniqueness_length_;
	}

	bool IsInnerShortEdge(set<VertexId> component, EdgeId edge) {
		return !IsUnique(edge) && component.count(graph_.EdgeStart(edge)) == 1
				&& component.count(graph_.EdgeEnd(edge)) == 1;
	}

	void ProcessShortEdge(FlowGraph<Graph> &fg, set<VertexId> component,
			EdgeId edge) {
		if (!IsInnerShortEdge(component, edge)) {
			fg.AddEdge(graph_.EdgeStart(edge), graph_.EdgeEnd(edge));
		}
	}

	void ProcessSource(FlowGraph<Graph> &fg, set<VertexId> component,
			EdgeId edge) {
		if (IsPlausible(edge) || IsUnique(edge)) {
			fg.AddSource(graph_.EdgeEnd(edge), 1);
		}
	}

	void ProcessSink(FlowGraph<Graph> &fg, set<VertexId> component,
			EdgeId edge) {
		if (IsPlausible(edge) || IsUnique(edge)) {
			fg.AddSink(graph_.EdgeStart(edge), 1);
		}
	}

	void ConstructFlowGraph(FlowGraph<Graph> &fg, set<VertexId> component) {
		for (auto it = component.begin(); it != component.end(); ++it) {
			fg.AddVertex(*it);
		}
		for (auto it = component.begin(); it != component.end(); ++it) {
			VertexId vertex = *it;
			auto outgoing = graph_.OutgoingEdges(vertex);
			for (auto it_edge = outgoing.begin(); it_edge != outgoing.end();
					++it_edge) {
				EdgeId edge = *it_edge;
				ProcessShortEdge(fg, component, edge);
				ProcessSink(fg, component, edge);
			}
			auto incoming = graph_.IncomingEdges(vertex);
			for (auto it_edge = incoming.begin(); it_edge != incoming.end();
					++it_edge) {
				EdgeId edge = *it_edge;
				ProcessSource(fg, component, edge);
			}
		}
	}

public:
	MaxFlowECRemover(Graph& graph, size_t max_length, size_t uniqueness_length,
			size_t plausibility_length, EdgeRemover<Graph>& edge_remover) :
			graph_(graph), max_length_(max_length), uniqueness_length_(
					uniqueness_length), plausibility_length_(
					plausibility_length), edge_remover_(edge_remover) {
		VERIFY(uniqueness_length >= plausibility_length);
		VERIFY(plausibility_length > max_length);
	}

	bool RemoveEdges() {
		LongEdgesExclusiveSplitter<Graph> splitter(graph_, uniqueness_length_);
		for (LongEdgesExclusiveSplitter<Graph> splitter(graph_,
				uniqueness_length_); splitter.Finished();) {
			auto component_vector = splitter.NextComponent();
			set<VertexId> component(component_vector.begin(),
					component_vector.end());
			FlowGraph<Graph> fg;
			ConstructFlowGraph(fg, component);
			MaxFlowFinder<Graph> mf_finder(fg);
			mf_finder.Find();
			if (!CheckCompleteFlow(fg)) {
				TRACE("Suspicious component! No edge delition!");
				continue;
			}
			StroglyConnectedComponentFinder<FlowGraph<Graph> > component_finder(
					fg);
			map<typename FlowGraph<Graph>::VertexId, size_t> colouring =
					component_finder.ColourComponents();
			set<EdgeId> to_remove = CollectUnusedEdges(component, fg,
					colouring);
			for (SmartSetIterator<Graph, EdgeId> it(graph_, to_remove.begin(), to_remove.end());
					!it.IsEnd(); ++it) {
				edge_remover_.DeleteEdge(*it);
			}
		}
		Compressor<Graph>(graph_).CompressAllVertices();
		Cleaner<Graph>(graph_).Clean();
		return false;
	}
private:
	DECL_LOGGER("MaxFlowECRemover")
	;
};
}

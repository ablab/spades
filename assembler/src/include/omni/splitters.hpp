#pragma once

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard_base.hpp"
#include "graph_component.hpp"
#include "dijkstra_tools/dijkstra_helper.hpp"
#include "component_filters.hpp"

namespace omnigraph {


template<typename Element>
class JSIterator {
public:

    virtual Element Next() = 0;

    virtual bool HasNext() = 0;

    virtual ~JSIterator() {
    }
};

template<class Graph>
class GraphSplitter : public JSIterator<GraphComponent<Graph>>{
private:
    const Graph& graph_;
public:
    GraphSplitter(const Graph& graph)
            : graph_(graph) {
    }

    const Graph& graph() const {
        return graph_;
    }
};

template<class Graph>
class PrecountedComponentSplitter : public GraphSplitter<Graph> {
    bool HasNext_;
    GraphComponent<Graph> component_;
public:

    template<class It>
    PrecountedComponentSplitter(const Graph &graph, It begin, It end)
            : GraphSplitter<Graph>(graph), HasNext_(false),
              component_(graph, begin, end) {
    }

    template<class It>
    PrecountedComponentSplitter(GraphComponent<Graph> component)
            : GraphSplitter<Graph>(component.g()), HasNext_(false),
              component_(component) {
    }

    GraphComponent<Graph> Next() {
        HasNext_ = false;
        return component_;
    }

//  virtual bool CheckPutVertex(VertexId /*vertex*/, EdgeId edge, size_t /*length*/) const {
//    return edges_.count(edge) != 0;
//  }
    bool HasNext() {
        return HasNext_;
    }
};

template<typename Element>
class RelaxingIterator : public JSIterator<Element> {
public:
    template<typename It>
    void Relax(It begin, It end) {
        Relax(vector<Element>(begin, end));
    }

//  virtual bool CheckProcessVertex(VertexId /*vertex*/, size_t distance) {
//    return distance <= bound_;
//  }
    virtual void Relax(const vector<Element> &v) = 0;

    virtual void Relax(Element) = 0;

    virtual ~RelaxingIterator() {
    }
};

template<class Collection>
class CollectionIterator : public RelaxingIterator<typename Collection::value_type> {
private:
    typedef typename Collection::value_type Element;
    typedef typename Collection::const_iterator Iter;
    shared_ptr<Collection> storage_;
    Iter current_;
    const Iter end_;
    set<Element> relaxed_;
public:
    CollectionIterator(const Collection &collection)
            : current_(collection.begin()), end_(collection.end()) {
    }

//  virtual bool CheckPutVertex(VertexId vertex, EdgeId /*edge*/, size_t /*length*/) const {
//    return subgraph_.count(vertex) != 0;
//  }
    CollectionIterator(shared_ptr<Collection> collection)
            : storage_(collection), current_(collection->begin()), end_(collection->end()) {
    }

    CollectionIterator(Iter begin, Iter end)
            : current_(begin), end_(end) {
    }

    Element Next() {
        if(!HasNext()) { //This function actually changes value of current! It is not just to verify!
            //fixme use VERIFY_MSG instead
            VERIFY(HasNext());
        }
        Element next = *current_;
        ++current_;
        return next;
    }

//public:
//  ErrorComponentSplitter(const Graph &graph, const set<EdgeId> &black_edges) :
//      base(graph), black_edges_(black_edges), iterator_(
//          graph.SmartEdgeBegin()) {
//    TRACE("ErrorComponentSplitter created and SmartIterator initialized");
//  }
//
//  virtual ~ErrorComponentSplitter() {
//  }
//
//  vector<VertexId> FindComponent(VertexId start_vertex) {
//    ComponentFinder<Graph> cf(this->graph(), black_edges_);
//    cf.run(start_vertex);
//    return cf.ReachedVertices();
//  }
//
//  vector<VertexId> FindNeighbourhood(VertexId start, size_t bound) {
//    NeighbourhoodFinder<Graph> nf(this->graph(), black_edges_, bound);
//    nf.run(start);
//    return nf.ReachedVertices();
//  }
//
//  size_t FindDiameter(const vector<VertexId> &component) {
//    set < VertexId > component_set(component.begin(), component.end());
//    size_t result = 0;
//    VertexId current = *(component.begin());
//    for (size_t i = 0; i < 4; i++) {
//      pair<VertexId, size_t> next = GetFarthest(current, component_set);
//      current = next.first;
//      result = next.second;
//    }
//    return result;
//  }
//
//  pair<VertexId, size_t> GetFarthest(VertexId v,
//      const set<VertexId> &component) {
//    SubgraphDijkstra<Graph> sd(this->graph(), component);
//    sd.run(v);
//    pair<VertexId, size_t> result(v, 0);
//    auto bounds = sd.GetDistances();
//    for (auto it = bounds.first; it != bounds.second; ++it) {
//      if (it->second > result.second) {
//        result = *it;
//      }
//    }
//    return result;
//  }
//
//  virtual vector<VertexId> NextComponent() {
//    TRACE("Construction of next component started");
//    if (Finished()) {
//      VERIFY(false);
//      return vector<VertexId>();
//    }
//    EdgeId next = *iterator_;
//    ++iterator_;
//    vector < VertexId > component = FindComponent(
//        this->graph().EdgeEnd(next));
//    TRACE("Error edges component constructed. It contains "
//            << component.size() << " vertices");
//    size_t component_size = FindDiameter(component);
//    TRACE("Diameter of component is " << component_size);
//    vector < VertexId > neighbourhood = FindNeighbourhood(
//        this->graph().EdgeEnd(next), (size_t) math::round(1.5 * (double) component_size));
//    TRACE("Error edges component neighborhood constructed. It contains "
//            << neighbourhood.size() << " vertices");
//    visited_.insert(component.begin(), component.end());
//    return neighbourhood;
//  }
//
//  virtual bool Finished() {
//    while (!iterator_.IsEnd()) {
//      if (black_edges_.find(*iterator_) != black_edges_.end()
//          && visited_.find(this->graph().EdgeEnd(*iterator_))
//              == visited_.end()) {
//        return false;
//      }
//      ++iterator_;
//    }
//    return true;
//  }
    bool HasNext() {
        while(current_ != end_ && relaxed_.count(*current_) == 1) {
            ++current_;
        }
        return current_ != end_;
    }

    void Relax(Element e) {
        relaxed_.insert(e);
    }

//template<class Graph>
//class ShortEdgeComponentNeighbourhoodFinder: public UnorientedDijkstra<Graph> {
//private:
//  typedef UnorientedDijkstra<Graph> base;
//protected:
//  typedef typename base::VertexId VertexId;
//  typedef typename base::EdgeId EdgeId;
//  typedef typename base::DistanceType distance_t;
//private:
//  distance_t bound_;
//public:
//  ShortEdgeComponentNeighbourhoodFinder(const Graph &graph, distance_t bound) :
//      UnorientedDijkstra<Graph>(graph), bound_(bound) {
//  }
//
//  virtual bool CheckProcessVertexVertexId (VertexId /*vertex*/, distance_t distance) {
//    return distance == 0;
//  }
//
//  virtual distance_t GetLength(EdgeId edge) const {
//    if (this->graph().length(edge) <= bound_)
//      return 0;
//    else
//      return 1;
//  }
    void Relax(const vector<Element> &v) {
        for (auto it = v.begin(); it != v.end(); ++it)
            Relax(*it);
    }

    virtual ~CollectionIterator() {
    }
};

template<class Graph>
class PathIterator : public RelaxingIterator<typename Graph::VertexId> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    vector<VertexId> path_;
    size_t current_;

    static vector<VertexId> ExtractVertices(const Graph &graph, const Path<EdgeId> &path) {
        vector<VertexId> result;
        for(size_t i = 0; i < path.size(); i++) {
        	if(i == 0 || path[i] != path[i - 1]) {
        		result.push_back(graph.EdgeStart(path[i]));
        		result.push_back(graph.EdgeEnd(path[i]));
        	}
        }
        return result;
    }

public:
    PathIterator(const Graph &graph, const Path<EdgeId> &path)
            : graph_(graph), path_(ExtractVertices(graph, path)), current_(0) {
    }

    VertexId Next() {
        if(!HasNext()) {
            VERIFY(HasNext());
        }
        VertexId next = path_[current_];
        Relax(next);
        return next;
    }

    bool HasNext() {
        return current_ < path_.size();
    }

    void Relax(const vector<VertexId> &v) {
        set<VertexId> toRelax(v.begin(), v.end());
        while(toRelax.count(path_[current_]) == 1)
            current_++;
    }

//public:
//  CountingDijkstra(const Graph &graph, size_t max_size,
//      size_t edge_length_bound) :
//      base(graph), max_size_(max_size), edge_length_bound_(
//          edge_length_bound), current_(0) {
//  }
//
//  virtual bool CheckPutVertex(VertexId /*vertex*/, EdgeId edge,
//      distance_t /*length*/) const {
//    if (current_ < max_size_) {
//      ++current_;
//    }
//    if (current_ < max_size_ && GetLength(edge) < inf) {
//      return true;
//    }
//    return false;
//  }
//
//  virtual bool CheckProcessVertex(VertexId /*vertex*/, distance_t /*distance*/) {
//    return current_ < max_size_;
//  }
//
//  virtual void init(VertexId /*start*/) {
//    current_ = 0;
//  }
//
//  virtual size_t GetLength(EdgeId edge) const {
//    if (this->graph().length(edge) <= edge_length_bound_)
//      //todo change back
////			return 1;
//      return this->graph().length(edge);
//    else
//      return inf;
//  }
    void Relax(VertexId e) {
        Relax(vector<VertexId>({e}));
    }
};

template<class Graph>
class AbstractNeighbourhoodFinder {
private:
    const Graph &graph_;
public:
    AbstractNeighbourhoodFinder(const Graph &graph) : graph_(graph) {
    }

    const Graph &graph() const {
        return graph_;
    }

    virtual GraphComponent<Graph> Find(typename Graph::VertexId v) = 0;

    virtual vector<typename Graph::VertexId> InnerVertices(const GraphComponent<Graph> &component) = 0;

    virtual ~AbstractNeighbourhoodFinder() {
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
    ComponentCloser(const Graph &graph, size_t edge_length_bound)
            : graph_(graph),
              edge_length_bound_(edge_length_bound) {
    }

    void CloseComponent(set<VertexId> &component) const {
        set<VertexId> additional_vertices;
        for (auto it = component.begin(); it != component.end(); ++it) {
            FOREACH (EdgeId e, graph_.OutgoingEdges(*it)) {
                if (graph_.length(e) >= edge_length_bound_) {
                    additional_vertices.insert(graph_.EdgeEnd(e));
                }
            }
            FOREACH (EdgeId e, graph_.IncomingEdges(*it)) {
                if (graph_.length(e) >= edge_length_bound_) {
                    additional_vertices.insert(graph_.EdgeStart(e));
                }
            }
        }
        component.insert(additional_vertices.begin(),
                         additional_vertices.end());
    }

    GraphComponent<Graph> CloseComponent(const GraphComponent<Graph>& component) const {
        set<VertexId> vertices(component.v_begin(), component.v_end());
        CloseComponent(vertices);
        return GraphComponent<Graph>(graph_, vertices.begin(), vertices.end());
    }
};

//This method finds a neighbourhood of a set of vertices. Vertices that are connected by an edge of length more than 600 are not considered as adjacent.
template<class Graph>
class ReliableNeighbourhoodFinder : public AbstractNeighbourhoodFinder<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    set<VertexId> FindNeighbours(const set<VertexId> &s) {
        set<VertexId> result(s.begin(), s.end());
        FOREACH (VertexId v, result) {
            FOREACH (EdgeId e, this->graph().AdjacentEdges(v)) {
                if(this->graph().length(e) <= edge_length_bound_) {
                    result.insert(this->graph().EdgeEnd(e));
                    result.insert(this->graph().EdgeStart(e));
                }
            }
        }
        return result;
    }

    set<VertexId> FindNeighbours(const set<VertexId> &s, size_t eps) {
        set<VertexId> result = s;
        for(size_t i = 0; i < eps; i++) {
            result = FindNeighbours(result);
        }
        return result;
    }

    set<VertexId> FindBorder(const GraphComponent<Graph> component) {
        set<VertexId> result;
        for(auto it = component.vertices().begin(); it != component.vertices().end(); ++it) {
            if(component.IsBorder(*it)) {
                result.insert(*it);
            }
        }
        return result;
    }

public:
    static const size_t DEFAULT_EDGE_LENGTH_BOUND = 500;
    static const size_t DEFAULT_MAX_SIZE = 100;

    const size_t edge_length_bound_;
    const size_t max_size_;

    ReliableNeighbourhoodFinder(const Graph &graph, size_t edge_length_bound =
                                        DEFAULT_EDGE_LENGTH_BOUND,
                                size_t max_size = DEFAULT_MAX_SIZE)
            : AbstractNeighbourhoodFinder<Graph>(graph),
              edge_length_bound_(edge_length_bound),
              max_size_(max_size) {
    }

//template<class Graph>
//class ReliableSplitterAlongPath: public GraphSplitter<Graph> {
//private:
//  typedef GraphSplitter<Graph> base;
//  typedef typename Graph::VertexId VertexId;
//  typedef typename Graph::EdgeId EdgeId;
//  size_t max_size_;
//  size_t edge_length_bound_;
//  set<VertexId> last_component_;
//  size_t current_index_;
//  MappingPath<EdgeId> path_;
//  Range covered_range_;
//  bool start_processed_;
//
//  //todo edge not used in the body
//  bool EdgeCovered(EdgeId /*edge*/) {
//    return last_component_.count(
//        this->graph().EdgeStart(path_[current_index_].first)) == 1
//        && last_component_.count(
//            this->graph().EdgeEnd(path_[current_index_].first)) == 1;
//  }
//
//  void SkipVisited() {
//    covered_range_.start_pos =
//        path_[current_index_].second.initial_range.start_pos;
//    covered_range_.end_pos =
//        path_[current_index_].second.initial_range.end_pos;
//    //always go forward at least one path element
//    ++current_index_;
//    while (current_index_ < path_.size()
//        && EdgeCovered(path_[current_index_].first)) {
//      covered_range_.end_pos =
//          path_[current_index_].second.initial_range.end_pos;
//      ++current_index_;
//    }
//  }
    GraphComponent<Graph> Find(typename Graph::VertexId v) {
    	auto cd = DijkstraHelper<Graph>::CreateCountingDijkstra(this->graph(), max_size_,
    			edge_length_bound_);
        cd.run(v);
        vector<VertexId> result_vector = cd.ReachedVertices();
        set<VertexId> result(result_vector.begin(), result_vector.end());
        ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
        cc.CloseComponent(result);
        return GraphComponent<Graph>(this->graph(), result.begin(),
                                     result.end());
    }

    vector<VertexId> InnerVertices(const GraphComponent<Graph> &component) {
        set<VertexId> border = FindNeighbours(FindBorder(component), 2);
        std::vector<VertexId> result;
        std::set_difference(component.vertices().begin(), component.vertices().end(), border.begin(), border.end(), std::inserter(result, result.end()));
        return vector<VertexId>(result.begin(), result.end());
    }
};

template<class Graph>
class ShortEdgeComponentFinder : public AbstractNeighbourhoodFinder<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    static const size_t DEFAULT_EDGE_LENGTH_BOUND = 100;

    const size_t edge_length_bound_;

    ShortEdgeComponentFinder(const Graph &graph, size_t edge_length_bound = DEFAULT_EDGE_LENGTH_BOUND)
            : AbstractNeighbourhoodFinder<Graph>(graph),
              edge_length_bound_(edge_length_bound) {
    }

    GraphComponent<Graph> Find(typename Graph::VertexId v) {
    	auto cd = DijkstraHelper<Graph>::CreateShortEdgeDijkstra(this->graph(), edge_length_bound_);
        cd.run(v);
        set<VertexId> result = cd.ProcessedVertices();
        return GraphComponent<Graph>(this->graph(), result.begin(),
                                     result.end());
    }

    vector<VertexId> InnerVertices(const GraphComponent<Graph> &component) {
        return vector<VertexId>(component.v_begin(), component.v_end());
    }
};

template<class Graph>
class FilteringSplitterWrapper : public GraphSplitter<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    shared_ptr<GraphComponentFilter<Graph>> checker_;
    boost::optional<GraphComponent<Graph>> next_;
public:
    FilteringSplitterWrapper(
            shared_ptr<GraphSplitter<Graph>> inner_splitter,
            shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()), inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
            VERIFY(false);
            return omnigraph::GraphComponent<Graph>(this->graph());
        }
        GraphComponent<Graph> result = next_.get();
        next_ = boost::optional<GraphComponent<Graph>>();
        return result;
    }

    bool HasNext() {
        while (!next_ && inner_splitter_->HasNext()) {
            GraphComponent<Graph> ne = inner_splitter_->Next();
            if (checker_->Check(ne)) {
                next_ = ne;
            }
        }
        return next_;
    }
private:
    DECL_LOGGER("FilteringSplitterWrapper");
};

//TODO  split combined component into several.
template<class Graph>
class CollectingSplitterWrapper : public GraphSplitter<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    shared_ptr<GraphComponentFilter<Graph>> checker_;
    boost::optional<GraphComponent<Graph>> next_;
    set<VertexId> filtered_;
public:
    CollectingSplitterWrapper(
            shared_ptr<GraphSplitter<Graph>> inner_splitter,
            shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()), inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
       		VERIFY(false);
           	return omnigraph::GraphComponent<Graph>(this->graph());
        } else {
        	if(next_) {
        		GraphComponent<Graph> result = next_.get();
        		next_ = boost::optional<GraphComponent<Graph>>();
        		return result;
        	} else {
           		GraphComponent<Graph> result(this->graph(), filtered_.begin(), filtered_.end(), false, "filtered");
           		filtered_.clear();
           		return result;
        	}
        }
    }

    bool HasNext() {
        while (!next_ && inner_splitter_->HasNext()) {
            GraphComponent<Graph> ne = inner_splitter_->Next();
            if (checker_->Check(ne)) {
                next_ = ne;
            } else {
            	filtered_.insert(ne.v_begin(), ne.v_end());
            }
        }
        return next_ || !filtered_.empty();
    }
private:
    DECL_LOGGER("FilteringSplitterWrapper");
};

template<class Graph>
class CondensingSplitterWrapper : public GraphSplitter<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    shared_ptr<GraphComponentFilter<Graph>> checker_;
    boost::optional<GraphComponent<Graph>> next_;

    string CutName(const string &name, size_t max_length) {
    	VERIFY(max_length >= 7);
    	size_t length = name.size();
    	if (length <= max_length)
    		return name;
    	else {
    		return name.substr(0, (max_length - 5) / 2) + "....." + name.substr(length - (max_length - 5) / 2, (max_length - 5) / 2);
    	}
    }

    GraphComponent<Graph> ConstructComponent() {
    	GraphComponent<Graph> next = inner_splitter_->Next();
    	if (checker_->Check(next)) {
    		return next;
    	}
    	set<VertexId> vertices(next.v_begin(), next.v_end());
    	string name = next.name();
    	for(size_t i = 0; i < 10 && inner_splitter_->HasNext(); i++) {
			next = inner_splitter_->Next();
			if (checker_->Check(next)) {
				next_ = next;
				break;
			} else {
				vertices.insert(next.v_begin(), next.v_end());
				name += ";";
				name += next.name();
			}
		}
		return GraphComponent<Graph>(this->graph(), vertices.begin(), vertices.end(), CutName(name, 60));
    }

public:
    CondensingSplitterWrapper(
            shared_ptr<GraphSplitter<Graph>> inner_splitter,
            shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()), inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
            VERIFY(false);
            return omnigraph::GraphComponent<Graph>(this->graph());
        }
        if(next_) {
        	GraphComponent<Graph> result = next_.get();
        	next_ = boost::optional<GraphComponent<Graph>>();
        	return result;
        } else {
        	return ConstructComponent();
        }
    }

    bool HasNext() {
    	if(next_)
    		return true;
    	if(!inner_splitter_->HasNext())
    		return false;
    	return true;
    }
private:
    DECL_LOGGER("FilteringSplitterWrapper");
};

template<class Graph>
class NeighbourhoodFindingSplitter : public GraphSplitter<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator_;
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> neighbourhood_finder_;

public:
    NeighbourhoodFindingSplitter(
            const Graph& graph,
            shared_ptr<RelaxingIterator<VertexId>> inner_iterator,
            shared_ptr<AbstractNeighbourhoodFinder<Graph>> neighbourhood_finder)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(inner_iterator),
              neighbourhood_finder_(neighbourhood_finder) {
    }

    NeighbourhoodFindingSplitter(
            const Graph& graph,
            shared_ptr<RelaxingIterator<VertexId>> inner_iterator)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(inner_iterator),
              neighbourhood_finder_(
                      make_shared<ReliableNeighbourhoodFinder<Graph>>(graph)) {
    }

    NeighbourhoodFindingSplitter(const Graph& graph)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(
                      make_shared<CollectionIterator<set<VertexId>>>(graph.begin(), graph.end())),
                      neighbourhood_finder_(make_shared<ReliableNeighbourhoodFinder<Graph>>(graph)) {
    }

    GraphComponent<Graph> Next() {
        VertexId next_vertex = inner_iterator_->Next();
        GraphComponent<Graph> result = neighbourhood_finder_->Find(next_vertex);
        vector<VertexId> to_relax = neighbourhood_finder_->InnerVertices(result);
        to_relax.push_back(next_vertex);
        inner_iterator_->Relax(to_relax);
        return result;
    }

    bool HasNext() {
        return inner_iterator_->HasNext();
    }
};

template<class Graph>
shared_ptr<GraphSplitter<Graph>> ReliableSplitter(const Graph &graph,
                            size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND,
                            size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE) {
    typedef typename Graph::VertexId VertexId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<CollectionIterator<set<VertexId>>>(graph.begin(), graph.end());
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph,
            inner_iterator, nf);
}

template<class Graph>
shared_ptr<GraphSplitter<Graph>> ReliableSplitterAlongPath(
        const Graph &graph, const Path<typename Graph::EdgeId>& path,
                                    size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND,
                                    size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE) {
    typedef typename Graph::VertexId VertexId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<
            PathIterator<Graph>>(graph, path);
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<
    		ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph,
                                                            inner_iterator, nf);
}

template<class Graph>
shared_ptr<GraphSplitter<Graph>> LongEdgesExclusiveSplitter(
        const Graph &graph, size_t bound =
                ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<
            CollectionIterator<set<VertexId>>>(graph.begin(), graph.end());
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<
            ShortEdgeComponentFinder<Graph>>(graph, bound);
    return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph,
                                                            inner_iterator, nf);
}

template<class Graph, typename Collection>
shared_ptr<GraphSplitter<Graph>> StandardSplitter(
        const Graph &graph, const Collection &collection, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<CollectionIterator<Collection>>(collection);
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<
            ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound,
                                                max_size);
    return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph, typename Collection>
shared_ptr<GraphSplitter<Graph>> StandardSplitter(
        const Graph &graph, shared_ptr<Collection> collection, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<CollectionIterator<Collection>>(collection);
    shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<
            ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound,
                                                max_size);
    return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph>
shared_ptr<GraphSplitter<Graph>> WholeGraphSplitter(
        const Graph &graph, size_t max_size,
        size_t edge_length_bound) {
    return NeighbourhoodFindingSplitter<Graph>(graph, graph.vertices(), max_size, edge_length_bound);
}

template<class Graph>
GraphComponent<Graph> VertexNeighborhood(
        const Graph &graph, typename Graph::VertexId vertex, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    vector<typename Graph::VertexId> vv = {vertex};
    shared_ptr<vector<typename Graph::VertexId>> sh_vv = make_shared<vector<typename Graph::VertexId>>(vv);
    return StandardSplitter<Graph>(graph, sh_vv, max_size, edge_length_bound)->Next();
}

//TODO make a method that draws a picture that contains given set of edges for sure. ? mb refactor this into just drawing instead of splitting?
template<class Graph>
GraphComponent<Graph> EdgeNeighborhood(
        const Graph &graph, typename Graph::EdgeId edge, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    vector<typename Graph::VertexId> vv = {graph.EdgeStart(edge)};
    shared_ptr<vector<typename Graph::VertexId>> sh_vv = make_shared<vector<typename Graph::VertexId>>(vv);
    return StandardSplitter<Graph>(graph, sh_vv, max_size, edge_length_bound)->Next();
}

}

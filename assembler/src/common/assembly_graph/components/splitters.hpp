#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_component.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
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
protected:
    //todo remove after returning to optional
    std::unique_ptr<GraphComponent<Graph>> MakeUniquePtr(GraphComponent<Graph>&& component) const {
        return std::unique_ptr<GraphComponent<Graph>>(new GraphComponent<Graph>(std::move(component)));
    }

    //todo remove after returning to optional
    GraphComponent<Graph> GetValueAndReset(std::unique_ptr<GraphComponent<Graph>>& component_ptr) const {
        VERIFY(component_ptr);
        auto answer = std::move(*component_ptr);
        component_ptr = nullptr;
        return answer;
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
        Relax(std::vector<Element>(begin, end));
    }

//  virtual bool CheckProcessVertex(VertexId /*vertex*/, size_t distance) {
//    return distance <= bound_;
//  }
    virtual void Relax(const std::vector<Element> &v) = 0;

    virtual void Relax(Element) = 0;

    virtual ~RelaxingIterator() {
    }
};

template<class It>
class RangeIterator : public RelaxingIterator<typename std::iterator_traits<It>::value_type> {
private:
    typedef typename std::iterator_traits<It>::value_type Element;
    It current_, end_;
    std::set<Element> relaxed_;
public:
    RangeIterator(It begin, It end)
            : current_(std::move(begin)), end_(std::move(end)) {
    }

    Element Next() {
        if (!HasNext()) { //This function actually changes value of current! It is not just to verify!
            //fixme use VERIFY_MSG instead
            VERIFY(HasNext());
        }
        Element next = *current_;
        ++current_;
        return next;
    }

    bool HasNext() {
        while (current_ != end_ && relaxed_.count(*current_) == 1) {
            ++current_;
        }
        return current_ != end_;
    }

    void Relax(Element e) { relaxed_.insert(e); }

    void Relax(const std::vector<Element> &v) {
        for (const auto& e : v)
            Relax(e);
    }

    virtual ~RangeIterator() {}
};

template<class Collection>
class CollectionIterator : public RangeIterator<typename Collection::const_iterator> {
private:
    std::shared_ptr<Collection> storage_;
public:
    CollectionIterator(const Collection &collection)
            : CollectionIterator::RangeIterator(collection.begin(), collection.end()) {
    }

    CollectionIterator(std::shared_ptr<Collection> collection)
            : CollectionIterator::RangeIterator(collection->begin(), collection->end()),
              storage_(collection) {}
};

template<class Graph>
class PathIterator : public RelaxingIterator<typename Graph::VertexId> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    typedef std::vector<VertexId> Path;
    Path path_;
    size_t current_;

    static Path ExtractVertices(const Graph &graph, const std::vector<EdgeId> &path) {
        Path result;
        for(size_t i = 0; i < path.size(); i++) {
            if(i == 0 || path[i] != path[i - 1]) {
                result.push_back(graph.EdgeStart(path[i]));
                result.push_back(graph.EdgeEnd(path[i]));
            }
        }
        return result;
    }

public:
    PathIterator(const Graph &graph, const std::vector<EdgeId> &path)
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

    void Relax(const std::vector<VertexId> &v) {
        std::set<VertexId> toRelax(v.begin(), v.end());
        while(toRelax.count(path_[current_]) == 1)
            current_++;
    }

    void Relax(VertexId e) {
        Relax(std::vector<VertexId>{e});
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

    virtual GraphComponent<Graph> Find(typename Graph::VertexId v) const = 0;

    virtual std::vector<typename Graph::VertexId> InnerVertices(const GraphComponent<Graph> &component) const = 0;

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

    void CloseComponent(std::set<VertexId> &component) const {
        std::set<VertexId> additional_vertices;
        for (auto it = component.begin(); it != component.end(); ++it) {
            for (EdgeId e : graph_.OutgoingEdges(*it)) {
                if (graph_.length(e) >= edge_length_bound_) {
                    additional_vertices.insert(graph_.EdgeEnd(e));
                }
            }
            for (EdgeId e : graph_.IncomingEdges(*it)) {
                if (graph_.length(e) >= edge_length_bound_) {
                    additional_vertices.insert(graph_.EdgeStart(e));
                }
            }
        }
        component.insert(additional_vertices.begin(),
                         additional_vertices.end());
    }

    GraphComponent<Graph> CloseComponent(const GraphComponent<Graph>& component) const {
        std::set<VertexId> vertices(component.v_begin(), component.v_end());
        CloseComponent(vertices);
        return GraphComponent<Graph>::FromVertices(graph_, vertices);
    }
};

template<class Graph>
class HighCoverageComponentFinder : public AbstractNeighbourhoodFinder<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    class CoverageBoundedDFS {
    private:
        const Graph &graph_;
        const double coverage_bound_;
        const size_t edge_limit_;
        mutable size_t edge_summary_length_;
        const size_t edge_summary_length_limit_;

        void Find(EdgeId edge, std::set<EdgeId> &result) const {
            if (edge_summary_length_ > edge_summary_length_limit_)
                return;

            if (result.size() > edge_limit_)
                return;

            if (math::ls(graph_.coverage(edge), coverage_bound_))
                return;

            if (result.count(edge) || result.count(graph_.conjugate(edge)))
                return;

            edge_summary_length_ += graph_.length(edge);
            result.insert(edge);
            result.insert(graph_.conjugate(edge));

            VertexId v = graph_.EdgeEnd(edge);
            for (auto e : graph_.IncidentEdges(v))
                Find(e, result);

            v = graph_.EdgeStart(edge);
            for (auto e : graph_.IncidentEdges(v))
                Find(e, result);
        }

        // FIXME: Get rid of recursion, it's ugly!
        void Fill(EdgeId edge, phmap::flat_hash_set<EdgeId> &processed) const {
            if (edge_summary_length_ >= edge_summary_length_limit_)
                return;

            if (processed.size() > edge_limit_)
                return;

            if (math::ls(graph_.coverage(edge), coverage_bound_))
                return;

            if (processed.count(edge) || processed.count(graph_.conjugate(edge)))
                return;

            edge_summary_length_ += graph_.length(edge);
            processed.insert(edge);
            processed.insert(graph_.conjugate(edge));

            VertexId v = graph_.EdgeEnd(edge);
            for (auto e : graph_.IncidentEdges(v))
                Fill(e, processed);

            v = graph_.EdgeStart(edge);
            for (auto e : graph_.IncidentEdges(v))
                Fill(e, processed);
        }

    public:
        CoverageBoundedDFS(const Graph &graph, double coverage_bound, size_t edge_summary_limit,
                           size_t edge_limit = 500)
                : graph_(graph),
                  coverage_bound_(coverage_bound),
                  edge_limit_(edge_limit),
                  edge_summary_length_(0),
                  edge_summary_length_limit_(edge_summary_limit) {
        }

        std::set<EdgeId> Find(VertexId v) const {
            edge_summary_length_ = 0;
            std::set<EdgeId> result;
            for (auto e : graph_.OutgoingEdges(v)) {
                Find(e, result);
            }
            for (auto e : graph_.IncomingEdges(v)) {
                Find(e, result);
            }
            return result;
        }

        size_t Fill(VertexId v) const {
            edge_summary_length_ = 0;
            phmap::flat_hash_set<EdgeId> processed;
            for (auto e : graph_.OutgoingEdges(v))
                Fill(e, processed);
            for (auto e : graph_.IncomingEdges(v)) {
                Fill(e, processed);
            }
            return edge_summary_length_;
        }

        size_t CumulativeEdgeLength() const {
            return edge_summary_length_;
        }
    };


    const double coverage_bound_;
    CoverageBoundedDFS dfs_helper;

public:
    HighCoverageComponentFinder(const Graph &graph,
                                double max_coverage,
                                size_t edge_sum_limit = std::numeric_limits<size_t>::max())
            : AbstractNeighbourhoodFinder<Graph>(graph),
              coverage_bound_(max_coverage),
              dfs_helper(graph, max_coverage, edge_sum_limit) {
    }

    GraphComponent<Graph> Find(VertexId v) const {
        std::set<EdgeId> result = dfs_helper.Find(v);
        return GraphComponent<Graph>::FromEdges(this->graph(), result, false);
    }

    size_t CumulativeEdgeLength(VertexId v) const {
        size_t result = dfs_helper.Fill(v);
        DEBUG("Total edge length for vertex " << v.int_id() << " is " << result);
        return result;
    }

    std::vector<VertexId> InnerVertices(const GraphComponent<Graph> &component) const {
        return std::vector<VertexId>(component.v_begin(), component.v_end());
    }
};


//This method finds a neighbourhood of a set of vertices. Vertices that are connected by an edge of length more than 600 are not considered as adjacent.
template<class Graph>
class ReliableNeighbourhoodFinder : public AbstractNeighbourhoodFinder<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    std::set<VertexId> FindNeighbours(const std::set<VertexId> &s) const {
        std::set<VertexId> result(s.begin(), s.end());
        for (VertexId v : result) {
            for (EdgeId e : this->graph().IncidentEdges(v)) {
                if(this->graph().length(e) <= edge_length_bound_) {
                    result.insert(this->graph().EdgeEnd(e));
                    result.insert(this->graph().EdgeStart(e));
                }
            }
        }
        return result;
    }

    std::set<VertexId> FindNeighbours(const std::set<VertexId> &s, size_t eps) const {
        std::set<VertexId> result = s;
        for(size_t i = 0; i < eps; i++) {
            result = FindNeighbours(result);
        }
        return result;
    }

    std::set<VertexId> FindBorder(const GraphComponent<Graph> &component) const {
        std::set<VertexId> result;
        utils::insert_all(result, component.entrances());
        utils::insert_all(result, component.exits());
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

    GraphComponent<Graph> Find(typename Graph::VertexId v) const {
        auto cd = DijkstraHelper<Graph>::CreateCountingDijkstra(this->graph(), max_size_,
                edge_length_bound_);
        cd.Run(v);
        auto result_vector = cd.ReachedVertices();
        std::set<VertexId> result(result_vector.begin(), result_vector.end());
        ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
        cc.CloseComponent(result);
        return GraphComponent<Graph>::FromVertices(this->graph(), result);
    }

    std::vector<VertexId> InnerVertices(const GraphComponent<Graph> &component) const {
        std::set<VertexId> border = FindNeighbours(FindBorder(component), 2);
        std::vector<VertexId> result;
        std::set_difference(component.vertices().begin(), component.vertices().end(),
                            border.begin(), border.end(),
                            std::inserter(result, result.end()));
        return result;
    }
};

template<class Graph>
class PathNeighbourhoodFinder : public AbstractNeighbourhoodFinder<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    VertexId OtherEnd(EdgeId e, VertexId v) const {
        if (this->graph().EdgeStart(e) == v)
            return this->graph().EdgeEnd(e);
        else
            return this->graph().EdgeStart(e);
    }

    bool Go(VertexId v, size_t curr_depth, std::set<VertexId> &grey, std::set<VertexId> &black) const {
        //allows single vertex to be visited many times with different depth values
        TRACE("Came to vertex " << this->graph().str(v) << " on depth " << curr_depth);
        if (curr_depth >= max_depth_) {
            TRACE("Too deep");
            return true;
        }
        if (grey.size() >= max_size_) {
            TRACE("Too many vertices");
            return false;
        }

        TRACE("Started processing of vertex " << this->graph().str(v));
        grey.insert(v);

        TRACE("Sorting incident edges");
        std::vector<EdgeId> incident_path, incident_non_path;
        for (EdgeId e : this->graph().IncidentEdges(v)) {
            if (path_edges_.count(e) != 0) {
                /*condition not to go backward*/
                if (this->graph().EdgeStart(e) == v) {
                    incident_path.push_back(e);
                }
            } else {
                incident_non_path.push_back(e);
            }
        }

        for (EdgeId e : incident_non_path) {
            if (this->graph().length(e) > edge_length_bound_) {
                TRACE("Edge " << this->graph().str(e) << " is too long");
                continue;
            }
            TRACE("Going along edge " << this->graph().str(e));
            if (!Go(OtherEnd(e, v), curr_depth + 1, grey, black))
                return false;
        }

        TRACE("End processing of vertex " << this->graph().str(v));
        black.insert(v);

        for (EdgeId e : incident_path) {
            if (grey.count(OtherEnd(e, v)) != 0)
                continue;
            TRACE("Going along next path edge " << this->graph().str(e));
            if (!Go(OtherEnd(e, v), 0, grey, black))
                return false;
        }

        return true;
    }

public:
    static const size_t DEFAULT_EDGE_LENGTH_BOUND = 500;
    static const size_t DEFAULT_MAX_DEPTH = 2;
    static const size_t DEFAULT_MAX_SIZE = 20;

    std::set<EdgeId> path_edges_;
    const size_t edge_length_bound_;
    const size_t max_size_;
    const size_t max_depth_;

    mutable std::set<VertexId> last_inner_;

    PathNeighbourhoodFinder(const Graph &graph, const std::vector<EdgeId> &path,
                            size_t edge_length_bound = DEFAULT_EDGE_LENGTH_BOUND,
                            size_t max_size = DEFAULT_MAX_SIZE,
                            size_t max_depth = DEFAULT_MAX_DEPTH)
            : AbstractNeighbourhoodFinder<Graph>(graph),
              path_edges_(path.begin(), path.end()),
              edge_length_bound_(edge_length_bound),
              max_size_(max_size),
              max_depth_(max_depth) {
    }


    GraphComponent<Graph> Find(VertexId v) const {
        TRACE("Starting from vertex " << this->graph().str(v));
        last_inner_.clear();
        std::set<VertexId> grey, black;
        Go(v, 0, grey, black);
        last_inner_ = black;
        last_inner_.insert(v);
        ComponentCloser<Graph>(this->graph(), 0).CloseComponent(grey);
        return GraphComponent<Graph>::FromVertices(this->graph(), grey);
    }

    std::vector<VertexId> InnerVertices(const GraphComponent<Graph> &/*component*/) const {
        return std::vector<VertexId>(last_inner_.begin(), last_inner_.end());
    }
private:
    DECL_LOGGER("PathNeighbourhoodFinder");
};

//todo delete and think if we really need hierarchy
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

    GraphComponent<Graph> Find(VertexId v) const {
        auto cd = DijkstraHelper<Graph>::CreateShortEdgeDijkstra(this->graph(), edge_length_bound_);
        cd.Run(v);
        return GraphComponent<Graph>::FromVertices(this->graph(), cd.ProcessedVertices());
    }

    std::vector<VertexId> InnerVertices(const GraphComponent<Graph> &component) const {
        return std::vector<VertexId>(component.v_begin(), component.v_end());
    }
};

template<class Graph>
class FilteringSplitterWrapper : public GraphSplitter<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    std::shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    std::shared_ptr<GraphComponentFilter<Graph>> checker_;
    std::shared_ptr<GraphComponent<Graph>> next_;
public:
    FilteringSplitterWrapper(
            std::shared_ptr<GraphSplitter<Graph>> inner_splitter,
            std::shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()),
              inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
            VERIFY(false);
            return omnigraph::GraphComponent<Graph>(this->graph());
        }
        auto result = next_;
        next_ = nullptr;
        return *result;
    }

    bool HasNext() {
        while (!next_ && inner_splitter_->HasNext()) {
            next_ = std::make_shared(inner_splitter_->Next());
            if (!checker_->Check(*next_)) {
                next_ = nullptr;
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

    std::shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    std::shared_ptr<GraphComponentFilter<Graph>> checker_;
    std::unique_ptr<GraphComponent<Graph>> next_;
    std::set<VertexId> filtered_;
public:
    CollectingSplitterWrapper(
            std::shared_ptr<GraphSplitter<Graph>> inner_splitter,
            std::shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()), inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
           VERIFY(false);
           return omnigraph::GraphComponent<Graph>::Empty(this->graph());
        } else {
            if (next_) {
                return this->GetValueAndReset(next_);
            } else {
                auto result = GraphComponent<Graph>::FromVertices(this->graph(),
                                                                  filtered_,
                                                                  false, "filtered");
                filtered_.clear();
                return result;
            }
        }
    }

    bool HasNext() {
        while (!next_ && inner_splitter_->HasNext()) {
            next_ = this->MakeUniquePtr(inner_splitter_->Next());
            if (!checker_->Check(*next_)) {
                filtered_.insert(next_->v_begin(), next_->v_end());
                next_ = nullptr;
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

    std::shared_ptr<GraphSplitter<Graph>> inner_splitter_;
    std::shared_ptr<GraphComponentFilter<Graph>> checker_;
    std::unique_ptr<GraphComponent<Graph>> next_;

    std::string CutName(const std::string &name, size_t max_length) {
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
        std::set<VertexId> vertices(next.v_begin(), next.v_end());
        std::string name = next.name();
        for(size_t i = 0; i < 10 && inner_splitter_->HasNext(); i++) {
            next = inner_splitter_->Next();
            if (checker_->Check(next)) {
                VERIFY(!next_);
                next_ = this->MakeUniquePtr(std::move(next));
                break;
            } else {
                vertices.insert(next.v_begin(), next.v_end());
                if (next.name() != "") {
                    name += ";";
                    name += next.name();
                }
            }
        }
        return GraphComponent<Graph>::FromVertices(this->graph(), vertices, false, CutName(name, 60));
    }


public:
    CondensingSplitterWrapper(
            std::shared_ptr<GraphSplitter<Graph>> inner_splitter,
            std::shared_ptr<GraphComponentFilter<Graph>> checker)
            : GraphSplitter<Graph>(inner_splitter->graph()), inner_splitter_(inner_splitter),
              checker_(checker) {
    }

    GraphComponent<Graph> Next() {
        if (!HasNext()) {
            VERIFY(false);
            return GraphComponent<Graph>(this->graph());
        }

        if (next_) {
            return this->GetValueAndReset(next_);
        } else {
            return ConstructComponent();
        }
    }

    bool HasNext() {
        if (next_)
            return true;
        if (!inner_splitter_->HasNext())
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
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator_;
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> neighbourhood_finder_;

public:
    NeighbourhoodFindingSplitter(
            const Graph& graph,
            std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator,
            std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> neighbourhood_finder)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(inner_iterator),
              neighbourhood_finder_(neighbourhood_finder) {
    }

    NeighbourhoodFindingSplitter(
            const Graph& graph,
            std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(inner_iterator),
              neighbourhood_finder_(
                      std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph)) {
    }

    NeighbourhoodFindingSplitter(const Graph& graph)
            : GraphSplitter<Graph>(graph),
              inner_iterator_(
                      std::make_shared<CollectionIterator<std::set<VertexId>>>(graph.begin(), graph.end())),
                      neighbourhood_finder_(std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph)) {
    }

    GraphComponent<Graph> Next() {
        VertexId next_vertex = inner_iterator_->Next();
        GraphComponent<Graph> result = neighbourhood_finder_->Find(next_vertex);
        auto to_relax = neighbourhood_finder_->InnerVertices(result);
        to_relax.push_back(next_vertex);
        inner_iterator_->Relax(to_relax);
        return result;
    }

    bool HasNext() {
        return inner_iterator_->HasNext();
    }
};

template<class Graph>
std::shared_ptr<GraphSplitter<Graph>> ReliableSplitter(const Graph &graph,
                            size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND,
                            size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator =
            std::make_shared<RangeIterator<typename Graph::VertexIt>>(graph.begin(), graph.end());
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf =
            std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph>
std::shared_ptr<GraphSplitter<Graph>> ConnectedSplitter(const Graph &graph,
                            size_t edge_length_bound = 1000000,
                            size_t max_size = 1000000) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator = std::make_shared<RangeIterator<typename Graph::VertexIt>>(graph.begin(), graph.end());
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph>
std::shared_ptr<GraphSplitter<Graph>> ReliableSplitterAlongPath(const Graph &graph,
                                                                const std::vector<typename Graph::EdgeId> &path,
                                                                size_t edge_length_bound = PathNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND,
                                                                size_t max_size = PathNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
                                                                size_t max_depth = PathNeighbourhoodFinder<Graph>::DEFAULT_MAX_DEPTH) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator = std::make_shared<PathIterator<Graph>>(graph, path);
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf =
            std::make_shared<PathNeighbourhoodFinder<Graph>>(graph, path, edge_length_bound, max_size, max_depth);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph>
std::shared_ptr<GraphSplitter<Graph>> LongEdgesExclusiveSplitter(const Graph &graph,
                                                                 size_t bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator =
            std::make_shared<RangeIterator<typename Graph::VertexIt>>(graph.begin(), graph.end());
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf =
            std::make_shared<ShortEdgeComponentFinder<Graph>>(graph, bound);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph, typename Collection>
std::shared_ptr<GraphSplitter<Graph>> StandardSplitter(
        const Graph &graph, const Collection &collection, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator =
            std::make_shared<CollectionIterator<Collection>>(collection);
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf =
            std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph, typename Collection>
std::shared_ptr<GraphSplitter<Graph>> StandardSplitter(
        const Graph &graph, std::shared_ptr<Collection> collection, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    typedef typename Graph::VertexId VertexId;
    std::shared_ptr<RelaxingIterator<VertexId>> inner_iterator =
            std::make_shared<CollectionIterator<Collection>>(collection);
    std::shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf =
            std::make_shared<ReliableNeighbourhoodFinder<Graph>>(graph, edge_length_bound, max_size);
    return std::make_shared<NeighbourhoodFindingSplitter<Graph>>(graph, inner_iterator, nf);
}

template<class Graph>
std::shared_ptr<GraphSplitter<Graph>> WholeGraphSplitter(const Graph &graph, size_t max_size, size_t edge_length_bound) {
    return NeighbourhoodFindingSplitter<Graph>(graph, graph.vertices(), max_size, edge_length_bound);
}

template<class Graph>
GraphComponent<Graph> VertexNeighborhood(
        const Graph &graph, typename Graph::VertexId vertex, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    std::vector<typename Graph::VertexId> vv = {vertex};
    auto sh_vv = std::make_shared<std::vector<typename Graph::VertexId>>(vv);
    return StandardSplitter<Graph>(graph, sh_vv, max_size, edge_length_bound)->Next();
}

//TODO make a method that draws a picture that contains given set of edges for sure. ? mb refactor this into just drawing instead of splitting?
template<class Graph>
GraphComponent<Graph> EdgeNeighborhood(
        const Graph &graph, typename Graph::EdgeId edge, size_t max_size = ReliableNeighbourhoodFinder<Graph>::DEFAULT_MAX_SIZE,
        size_t edge_length_bound = ReliableNeighbourhoodFinder<Graph>::DEFAULT_EDGE_LENGTH_BOUND) {
    std::vector<typename Graph::VertexId> vv = {graph.EdgeStart(edge)};
    auto sh_vv = std::make_shared<std::vector<typename Graph::VertexId>>(vv);
    return StandardSplitter<Graph>(graph, sh_vv, max_size, edge_length_bound)->Next();
}

}

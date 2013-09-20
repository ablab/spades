#pragma once

#include "standard_base.hpp"
#include "graph_component.hpp"
#include "graph_processing_algorithm.hpp"

namespace omnigraph {

namespace simplification {

template<class EdgeContainer>
void SingleEdgeAdapter(
        const EdgeContainer& edges,
        boost::function<void(typename EdgeContainer::value_type)> single_edge_handler_f) {
    FOREACH(auto e, edges) {
        single_edge_handler_f(e);
    }
}

template<class Graph>
void VisualizeNontrivialComponentAutoInc(
        const Graph& g, const set<typename Graph::EdgeId>& edges,
        const string& folder, const GraphLabeler<Graph>& labeler,
        shared_ptr<visualization::GraphColorer<Graph>> colorer) {
    static size_t cnt = 0;
    if (edges.size() > 1) {
        set<typename Graph::VertexId> vertices;
        FOREACH(auto e, edges) {
            vertices.insert(g.EdgeStart(e));
            vertices.insert(g.EdgeEnd(e));
        }
        visualization::WriteComponent(
                GraphComponent<Graph>(g, vertices.begin(), vertices.end()),
                folder + ToString(cnt++) + ".dot", colorer, labeler);
    }
}

namespace relative_coverage {

template<class Graph>
vector<typename Graph::EdgeId> AdjacentEdges(const Graph& g, typename Graph::VertexId v) {
    vector<typename Graph::EdgeId> answer;
    push_back_all(answer, g.OutgoingEdges(v));
    push_back_all(answer, g.IncomingEdges(v));
    return answer;
}

template<class Graph>
class Component {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    set<EdgeId> edges_;
    set<VertexId> inner_vertices_;
    set<VertexId> border_;
    set<VertexId> terminating_vertices_;
    //maybe use something more sophisticated in future
    size_t cumm_length_;
    bool contains_deadends_;

    //if edge start = edge end = v returns v
    VertexId OppositeEnd(EdgeId e, VertexId v) const {
        VERIFY(g_.EdgeStart(e) == v
                || g_.EdgeEnd(e) == v);
//      VERIFY(remover_.g.EdgeStart(e) != remover_.g.EdgeEnd(e));
        if (g_.EdgeStart(e) == v) {
            return g_.EdgeEnd(e);
        } else {
            return g_.EdgeStart(e);
        }
    }

    void RemoveFromBorder(VertexId v) {
        size_t cnt = border_.erase(v);
        VERIFY(cnt);
    }

public:

    Component(const Graph& g, EdgeId e) : g_(g), cumm_length_(0), contains_deadends_(false) {
        edges_.insert(e);
        cumm_length_ += g_.length(e);
        border_.insert(g.EdgeStart(e));
        border_.insert(g.EdgeEnd(e));
    }

    void MakeInner(VertexId v) {
//        INFO("Checking if vertex " << g_.str(v) << " is tip.");
        VERIFY(border_.count(v) > 0);
        if (g_.IsDeadEnd(v) || g_.IsDeadStart(v)) {
//            INFO("Tip, flag put");
            contains_deadends_ = true;
        }
        inner_vertices_.insert(v);
        FOREACH(EdgeId e, AdjacentEdges(g_, v)) {
            //seems to correctly handle loops
            if (edges_.count(e) == 0) {
                edges_.insert(e);
                cumm_length_ += g_.length(e);
                VertexId other_end = OppositeEnd(e, v);
                if (inner_vertices_.count(other_end) == 0) {
                    border_.insert(other_end);
                }
            }
        }
        RemoveFromBorder(v);
    }

    void TerminateOnVertex(VertexId v) {
        terminating_vertices_.insert(v);
        RemoveFromBorder(v);
    }

    VertexId NextBorderVertex() const {
        return *border_.begin();
    }

    bool IsBorderEmpty() const {
        return border_.empty();
    }

    const set<EdgeId>& edges() const {
        return edges_;
    }

    bool contains(EdgeId e) const {
        return edges_.count(e) > 0;
    }

    const set<VertexId>& terminating_vertices() const {
        return terminating_vertices_;
    }

    const Graph& g() const {
        return g_;
    }

    size_t inner_vertex_cnt() const {
        return inner_vertices_.size();
    }

    size_t length() const {
        return cumm_length_;
    }

    bool contains_deadends() const {
        return contains_deadends_;
    }
};

template<class Graph>
class ComponentChecker;

template<class Graph>
class ComponentSearcher;

//currently works with conjugate graphs only (due to the assumption in the outer cycle)
template<class Graph>
class RelativeCoverageComponentRemover : public EdgeProcessingAlgorithm<Graph> {
public:
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef boost::function<double(EdgeId, VertexId)> LocalCoverageFT;
    typedef typename ComponentRemover<Graph>::HandlerF HandlerF;
private:

    LocalCoverageFT local_coverage_f_;
    size_t length_bound_;
    double min_coverage_gap_;
    size_t tip_allowing_length_bound_;
    size_t longest_connecting_path_bound_;
    //todo use it as guarding threshold?
    double max_coverage_;
    //bound on the number of inner vertices
    size_t vertex_count_limit_;
    ComponentRemover<Graph> component_remover_;

    boost::function<bool(EdgeId)> edge_classifier_;

//  double DetailLocalCoverage(EdgeId e, VertexId v) const {
//    INFO("Local coverage of edge " << this->g().str(e) << " around vertex "
//        << this->g().str(v) << " was " << local_coverage_f_(e, v));
//    return local_coverage_f_(e, v);
//  }
//
//  double MaxLocalCoverage(EdgeId e, VertexId v) const {
//      return std::max(DetailLocalCoverage(e, v), this->graph().coverage(e));
//  }
//
//
//  double MinLocalCoverage(EdgeId e, VertexId v) const {
//      return std::min(DetailLocalCoverage(e, v), this->graph().coverage(e));
//  }

    double LocalCoverage(EdgeId e, VertexId v) const {
        INFO("Local coverage of edge " << this->g().str(e) << " around vertex " << this->g().str(v) << " was " << local_coverage_f_(e, v));
        return local_coverage_f_(e, v);
    }

    double MaxLocalCoverage(const vector<EdgeId>& edges, VertexId v) const {
        double answer = 0.0;
        FOREACH(EdgeId e, edges) {
            answer = max(answer, LocalCoverage(e, v));
        }
        return answer;
    }

    bool CheckAnyHighlyCovered(const vector<EdgeId>& edges, VertexId v,
                               double base_coverage) const {
        return math::gr(MaxLocalCoverage(edges, v),
                        base_coverage * min_coverage_gap_);
    }

    friend class ComponentSearcher<Graph>;

    double RelativeCoverageToReport(VertexId v, double base_coverage) const {
        return std::min(MaxLocalCoverage(this->g().OutgoingEdges(v), v),
                        MaxLocalCoverage(this->g().IncomingEdges(v), v))
                / base_coverage;
    }

public:
//todo make some useful order and stop condition
    RelativeCoverageComponentRemover(
            Graph& g, LocalCoverageFT local_coverage_f, size_t length_bound,
            double min_coverage_gap, size_t tip_allowing_length_bound,
            size_t longest_connecting_path_bound,
            double max_coverage = std::numeric_limits<size_t>::max(),
            HandlerF handler_function = 0, size_t vertex_count_limit = 10,
            boost::function<bool(EdgeId)> edge_classifier = 0)
            : base(g),
              local_coverage_f_(local_coverage_f),
              length_bound_(length_bound),
              min_coverage_gap_(min_coverage_gap),
              tip_allowing_length_bound_(tip_allowing_length_bound),
              longest_connecting_path_bound_(longest_connecting_path_bound),
              max_coverage_(max_coverage),
              vertex_count_limit_(vertex_count_limit),
              component_remover_(g, handler_function),
              edge_classifier_(edge_classifier) {
        VERIFY(math::gr(min_coverage_gap, 1.));
        VERIFY(tip_allowing_length_bound <= length_bound);
        INFO("Coverage gap " << min_coverage_gap_);
    }

    //todo change qualifiers
protected:

    /*virtual*/
    bool ProcessEdge(EdgeId e) {
        INFO("Processing edge " << this->g().str(e));

        //here we use that the graph is conjugate!
        VertexId v = this->g().EdgeStart(e);

        if (this->g().IsDeadEnd(v) || this->g().IsDeadStart(v)) {
            INFO("Tip or isolated");
            return false;
        }

        double local_cov = LocalCoverage(e, v);

        INFO("Local coverage around start " << this->g().str(v) << " is " << local_cov);

        //temporary
        if (edge_classifier_ && edge_classifier_(e)) {
            VertexId v2 = this->g().EdgeEnd(e);
            INFO("Chimeric edge. Relative coverage info: "
                    << std::min(RelativeCoverageToReport(v, LocalCoverage(e, v)), RelativeCoverageToReport(v2, LocalCoverage(e, v2)))
                    << " "
                    << std::max(RelativeCoverageToReport(v, LocalCoverage(e, v)), RelativeCoverageToReport(v2, LocalCoverage(e, v2))));
        }

        //since min_coverage_gap_ > 1, we don't need to think about e here
        INFO("Checking presence of highly covered edges around start")
        if (CheckAnyHighlyCovered(this->g().OutgoingEdges(v), v, local_cov)
                && CheckAnyHighlyCovered(this->g().IncomingEdges(v), v,
                                         local_cov)) {
            INFO("Looking for component");
            ComponentChecker<Graph> checker(this->g(), vertex_count_limit_, length_bound_,
                                            tip_allowing_length_bound_,
                                            longest_connecting_path_bound_);
            //case of e being loop is handled implicitly!
            ComponentSearcher<Graph> component_searcher(
                    *this, checker, e);
            if (component_searcher.FindComponent()) {
                INFO("Deleting component");
                const Component<Graph>& component = component_searcher.component();
                component_remover_.DeleteComponent(component.edges());
                return true;
            } else {
                INFO("Failed to find component");
            }
        } else {
            INFO("No highly covered edges around");
        }
        return false;
    }
private:
    DECL_LOGGER("RelativeCoverageComponentRemover")
    ;
};

template<class Graph>
class LongestPathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Component<Graph>& component_;
    const Graph& g_;
    map<VertexId, size_t> max_distance_;
    vector<VertexId> vertex_stack_;
    bool cycle_detected_;

    //-1u if can't be counted yet
    size_t TryGetMaxDistance(VertexId v) {
        if (max_distance_.count(v) > 0)
            return max_distance_[v];

        size_t answer = 0;
        FOREACH (EdgeId e, g_.IncomingEdges(v)) {
            VertexId start = g_.EdgeStart(e);
            if (component_.contains(e)) {
                if (max_distance_.count(start) == 0) {
                    if (std::find(vertex_stack_.begin(), vertex_stack_.end(), start) != vertex_stack_.end()) {
                        cycle_detected_ = true;
                    }
                    vertex_stack_.push_back(start);
                    return -1u;
                } else {
                    answer = std::max(answer, max_distance_[start] + g_.length(e));
                }
            }
        }
        return answer;
    }

    void ProcessVertex(VertexId init_v) {
        vertex_stack_.push_back(init_v);
        while (!vertex_stack_.empty()) {
            if (cycle_detected_)
                return;

            VertexId v = vertex_stack_.back();
            size_t max_dist = TryGetMaxDistance(v);
            if (max_dist != -1u) {
                max_distance_[v] = max_dist;
                vertex_stack_.pop_back();
            }
        }
    }

public:
    LongestPathFinder(const Component<Graph>& component)
    : component_(component), g_(component.g()), cycle_detected_(false) {
    }

    //-1u if component contains a cycle
    size_t Find() {
        size_t answer = 0;
        FOREACH(VertexId v, component_.terminating_vertices()) {
            ProcessVertex(v);
            if (cycle_detected_)
                return -1u;
            VERIFY(max_distance_.count(v) > 0);
            answer = std::max(answer, get(max_distance_, v));
        }
        return answer;
    }
};

template<class Graph>
class ComponentChecker {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    size_t vertex_count_limit_;
    size_t length_bound_;
    size_t tip_allowing_length_bound_;
    size_t longest_connecting_path_bound_;

public:
    ComponentChecker(const Graph& g, size_t vertex_count_limit, size_t length_bound,
                     size_t tip_allowing_length_bound,
                     size_t longest_connecting_path_bound)
            : g_(g), vertex_count_limit_(vertex_count_limit),
              length_bound_(length_bound),
              tip_allowing_length_bound_(tip_allowing_length_bound),
              longest_connecting_path_bound_(longest_connecting_path_bound) {
    }

    bool SizeCheck(const Component<Graph>& component) const {
        if (component.inner_vertex_cnt() > vertex_count_limit_) {
            INFO("Too many vertices! More than " << vertex_count_limit_);
            return false;
        }
        if (!component.contains_deadends()
                && component.length() > length_bound_) {
            INFO("Too long component of length " << component.length() << "! Longer than length bound " << length_bound_);
            return false;
        }
        if (component.length() > tip_allowing_length_bound_) {
            INFO("Too long component of length " << component.length() << "! Longer than tip allowing length bound " << tip_allowing_length_bound_);
            return false;
        }
        return true;
    }

    bool FullCheck(const Component<Graph>& component) const {
        size_t longest_connecting_path = LongestPathFinder<Graph>(component).Find();
        return SizeCheck(component) && (longest_connecting_path == -1u || longest_connecting_path < longest_connecting_path_bound_);
    }

private:
    DECL_LOGGER("RelativelyLowCoveredComponentChecker")
    ;
};

template<class Graph>
class ComponentSearcher {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const RelativeCoverageComponentRemover<Graph>& remover_;
    const ComponentChecker<Graph>& checker_;
    const Graph& g_;
    Component<Graph> component_;

public:
    ComponentSearcher(
            const RelativeCoverageComponentRemover<Graph>& remover,
            const ComponentChecker<Graph>& checker,
            EdgeId first_edge)
            : remover_(remover), checker_(checker),
              g_(remover_.g()), component_(g_, first_edge) {
    }

    bool FindComponent() {
        while (!component_.IsBorderEmpty() && checker_.SizeCheck(component_)) {
            VertexId v = component_.NextBorderVertex();

            INFO("Checking if vertex " << g_.str(v) << " is terminating.");
            //checking if there is a sufficient coverage gap
            if (!IsTerminateVertex(v)) {
                INFO("Not terminating, adding neighbourhood");
                component_.MakeInner(v);
            } else {
                INFO("Terminating");
                component_.TerminateOnVertex(v);
            }
        }
        return checker_.FullCheck(component_);
    }

    const Component<Graph>& component() const {
        return component_;
    }

private:

    bool IsTerminateVertex(VertexId v) const {
        double base_coverage = remover_.MaxLocalCoverage(
                RetainEdgesFromComponent(AdjacentEdges(g_, v)), v);
        return CheckAnyFilteredHighlyCovered(g_.OutgoingEdges(v),
                                             v, base_coverage)
                && CheckAnyFilteredHighlyCovered(
                        g_.IncomingEdges(v), v, base_coverage);
    }

    bool CheckAnyFilteredHighlyCovered(const vector<EdgeId>& edges,
                                       VertexId v,
                                       double base_coverage) const {
        return remover_.CheckAnyHighlyCovered(
                FilterEdgesFromComponent(edges), v, base_coverage);
    }

    vector<EdgeId> FilterEdgesFromComponent(
            const vector<EdgeId>& edges) const {
        vector<EdgeId> answer;
        FOREACH(EdgeId e, edges) {
            if (!component_.contains(e)) {
                answer.push_back(e);
            }
        }
        return answer;
    }

    vector<EdgeId> RetainEdgesFromComponent(
            const vector<EdgeId>& edges) const {
        vector<EdgeId> answer;
        FOREACH(EdgeId e, edges) {
            if (component_.contains(e)) {
                answer.push_back(e);
            }
        }
        return answer;
    }

    DECL_LOGGER("RelativelyLowCoveredComponentSearcher")
    ;
};

}
}

}

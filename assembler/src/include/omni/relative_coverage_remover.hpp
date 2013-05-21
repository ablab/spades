#pragma once

#include "standard_base.hpp"
#include "graph_component.hpp"
#include "graph_processing_algorithm.hpp"

namespace omnigraph {

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
    const GraphColorer<Graph>& colorer) {
  static size_t cnt = 0;
  if (edges.size() > 1) {
    set<typename Graph::VertexId> vertices;
    FOREACH(auto e, edges) {
      vertices.insert(g.EdgeStart(e));
      vertices.insert(g.EdgeEnd(e));
    }
    WriteComponent(GraphComponent<Graph>(g, vertices.begin(), vertices.end()),
                   folder + ToString(cnt++) + ".dot", colorer, labeler);
  }
}

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
  size_t max_coverage_;
  //bound on the number of inner vertices
  size_t vertex_count_limit_;
  ComponentRemover<Graph> component_remover_;

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
    INFO("Local coverage of edge " << this->g().str(e) << " around vertex "
        << this->g().str(v) << " was " << local_coverage_f_(e, v));
    return local_coverage_f_(e, v);
  }

  bool CheckAnyHighlyCovered(const vector<EdgeId>& edges, VertexId v,
                             double base_coverage) const {
    FOREACH(EdgeId e, edges) {
      //should be gr and not ge to deal with case of 0 local coverage in tests
      if (math::gr(LocalCoverage(e, v), base_coverage * min_coverage_gap_))
        return true;
    }
    return false;
  }

  class RelativelyLowCoveredComponentSearcher;

  friend class RelativelyLowCoveredComponentSearcher;

  class RelativelyLowCoveredComponentSearcher {
    const RelativeCoverageComponentRemover& remover_;

    set<EdgeId> component_;
    set<VertexId> inner_vertices_;
    set<VertexId> border_;

    //maybe use something more sophisticated in future
    size_t component_length_;
   public:
    RelativelyLowCoveredComponentSearcher(
        const RelativeCoverageComponentRemover& remover, EdgeId first_edge,
        VertexId first_border_vertex)
        : remover_(remover), component_length_(0) {
      component_.insert(first_edge);
      component_length_ += remover_.g().length(first_edge);
      border_.insert(first_border_vertex);
    }

    bool FindComponent() {
      while (!border_.empty()) {
        VertexId v = *border_.begin();
        border_.erase(v);

        //all the tips were already removed! It is dangerous to support them here.
        if (remover_.g().IsDeadEnd(v) || remover_.g().IsDeadStart(v))
          return false;

        INFO("Checking if vertex " << remover_.g().str(v) << " is terminating.");
        //checking if there is a sufficient coverage gap
        if (!IsTerminateVertex(v)) {
          INFO("Not terminating, adding neighbourhood");
          inner_vertices_.insert(v);
          FOREACH(EdgeId e, AdjacentEdges(v)) {
            //seems to correctly handle loops
            component_.insert(e);
            component_length_ += remover_.g().length(e);
            VertexId other_end = OppositeEnd(e, v);
            if (inner_vertices_.count(other_end) == 0) {
              border_.insert(other_end);
            }
          }
        } else {
          INFO("Terminating");
          //do nothing, we already erased v from the border
        }
        if (inner_vertices_.size() > remover_.vertex_count_limit_) {
          INFO("Too many vertices! More than " << remover_.vertex_count_limit_);
          return false;
        }
        if (component_length_ > remover_.length_bound_) {
          INFO("Too long component! Longer than " << remover_.length_bound_);
          return false;
        }
      }
      return true;
    }

    const set<EdgeId>& component() const {
      return component_;
    }

   private:

    bool IsTerminateVertex(VertexId v) const {
      double base_coverage = MaxLocalCoverage(
          RetainEdgesFromComponent(AdjacentEdges(v)), v);
      return CheckAnyFilteredHighlyCovered(remover_.g().OutgoingEdges(v), v,
                                           base_coverage)
          && CheckAnyFilteredHighlyCovered(remover_.g().IncomingEdges(v), v,
                                           base_coverage);
    }

    bool CheckAnyFilteredHighlyCovered(const vector<EdgeId>& edges, VertexId v,
                                       double base_coverage) const {
      return remover_.CheckAnyHighlyCovered(FilterEdgesFromComponent(edges), v,
                                            base_coverage);
    }

    //if edge start = edge end = v returns v
    VertexId OppositeEnd(EdgeId e, VertexId v) const {
      VERIFY(remover_.g().EdgeStart(e) == v || remover_.g().EdgeEnd(e) == v);
//      VERIFY(remover_.g.EdgeStart(e) != remover_.g.EdgeEnd(e));
      if (remover_.g().EdgeStart(e) == v) {
        return remover_.g().EdgeEnd(e);
      } else {
        return remover_.g().EdgeStart(e);
      }
    }

    vector<EdgeId> AdjacentEdges(VertexId v) const {
      vector<EdgeId> answer;
      push_back_all(answer, remover_.g().OutgoingEdges(v));
      push_back_all(answer, remover_.g().IncomingEdges(v));
      return answer;
    }

    vector<EdgeId> FilterEdgesFromComponent(const vector<EdgeId>& edges) const {
      vector<EdgeId> answer;
      FOREACH(EdgeId e, edges) {
        if (component_.count(e) == 0) {
          answer.push_back(e);
        }
      }
      return answer;
    }

    vector<EdgeId> RetainEdgesFromComponent(const vector<EdgeId>& edges) const {
      vector<EdgeId> answer;
      FOREACH(EdgeId e, edges) {
        if (component_.count(e) > 0) {
          answer.push_back(e);
        }
      }
      return answer;
    }

    double MaxLocalCoverage(const vector<EdgeId>& edges, VertexId v) const {
      double answer = 0.0;
      FOREACH(EdgeId e, edges) {
        answer = max(answer, remover_.LocalCoverage(e, v));
      }
      return answer;
    }

    DECL_LOGGER("RelativelyLowCoveredComponentSearcher");
  };

 public:
//todo make some useful order and stop condition
  RelativeCoverageComponentRemover(
      Graph& g, LocalCoverageFT local_coverage_f, size_t length_bound,
      double min_coverage_gap, size_t max_coverage =
          std::numeric_limits<size_t>::max(),
      HandlerF handler_function = 0, size_t vertex_count_limit = 10)
      : base(g),
        local_coverage_f_(local_coverage_f),
        length_bound_(length_bound),
        min_coverage_gap_(min_coverage_gap),
        max_coverage_(max_coverage),
        vertex_count_limit_(vertex_count_limit),
        component_remover_(g, handler_function) {
    VERIFY(math::gr(min_coverage_gap, 1.));
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
    //since min_coverage_gap_ > 1, we don't need to think about e here
    INFO("Checking presence of highly covered edges around start")
    if (CheckAnyHighlyCovered(this->g().OutgoingEdges(v), v, local_cov)
        && CheckAnyHighlyCovered(this->g().IncomingEdges(v), v, local_cov)) {
      INFO("Looking for component");
      //case of e being loop is handled implicitly!
      RelativelyLowCoveredComponentSearcher/*<Graph>*/component_searcher(
          *this, e, this->g().EdgeEnd(e));
      if (component_searcher.FindComponent()) {
        INFO("Deleting component");
        auto component = component_searcher.component();
        component_remover_.DeleteComponent(component);
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

}

#pragma once

#include "standard_base.hpp"
#include "graph_component.hpp"
#include "graph_processing_algorithm.hpp"

namespace omnigraph {

template<class Graph>
class RelativeCoverageComponentRemover : public EdgeProcessingAlgorithm<Graph> {
  typedef EdgeProcessingAlgorithm<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef boost::function<double(EdgeId, VertexId)> local_coverage_f_t;

  double min_coverage_gap_;
  local_coverage_f_t local_coverage_f_;
  size_t vertex_count_limit_;
  size_t max_coverage_;
  ComponentRemover<Graph> component_remover_;

  class RelativelyLowCoveredComponentSearcher;

  friend class RelativelyLowCoveredComponentSearcher;

  class RelativelyLowCoveredComponentSearcher {
    const RelativeCoverageComponentRemover& remover_;

    set<EdgeId> component_;
    set<VertexId> inner_vertices_;
    set<VertexId> border_;
   public:
    RelativelyLowCoveredComponentSearcher(
        const RelativeCoverageComponentRemover& remover)
        : remover_(remover) {

    }

    bool FindComponent() {
      while (!border_.empty()) {
        VertexId v = *border_.begin();
        //all the tips vere already removed
        if (remover_.g.IsDeadEnd(v) || remover_.g.IsDeadStart(v))
          return false;
        //checking if there is a sufficient coverage gap
        if (!IsTerminateVertex(v)) {
          border_.erase(v);
          inner_vertices_.insert(v);
          FOREACH(EdgeId e, AdjacentEdges(v)) {

          }
        }


      }
    }

   private:

    double LocalCoverage(EdgeId e, VertexId v) const {
      return remover_.local_coverage_f_(e, v);
    }

    bool IsTerminateVertex(VertexId v) const {
      double base_coverage = MaxLocalCoverage(
          RetainEdgesFromComponent(AdjacentEdges(v)), v);
      return CheckAnyFilteredHighlyCovered(remover_.g.OutgoingEdges(v), v,
                                           base_coverage)
          && CheckAnyFilteredHighlyCovered(remover_.g.IncomingEdges(v), v,
                                           base_coverage);
    }

    bool CheckAnyFilteredHighlyCovered(const vector<EdgeId>& edges, VertexId v,
                                       double base_coverage) {
      return CheckAnyHighlyCovered(FilterEdgesFromComponent(edges), v,
                                   base_coverage);
    }

    VertexId OppositeEnd(EdgeId e, VertexId v) const {
      VERIFY(remover_.g.EdgeStart(e) != remover_.g.EdgeEnd(e));
      if (remover_.g.EdgeStart(e) == v) {
        return remover_.g.EdgeEnd(e);
      } else if (remover_.g.EdgeEnd(e) == v) {
        return remover_.g.EdgeStart(e);
      } else {
        VERIFY(false);
        return 0.0;
      }
    }

    vector<EdgeId> AdjacentEdges(VertexId v) const {
      vector<EdgeId> answer;
      push_back_all(answer, remover_.g.OutgoingEdges(v));
      push_back_all(answer, remover_.g.IncomingEdges(v));
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

    bool CheckAnyHighlyCovered(const vector<EdgeId>& edges, VertexId v,
                               double base_coverage) const {
      FOREACH(EdgeId e, edges) {
        if (math::ge(LocalCoverage(e, v),
                     base_coverage * remover_.min_coverage_gap_))
          return true;
      }
      return false;
    }

    double MaxLocalCoverage(const vector<EdgeId>& edges, VertexId v) {
      double answer = 0.0;
      FOREACH(EdgeId e, edges) {
        answer = max(answer, LocalCoverage(e, v));
      }
      return answer;
    }

  };

 public:
//todo make some useful order and stop condition
  RelativeCoverageComponentRemover(
      Graph& g, double min_coverage_gap, local_coverage_f_t local_coverage_f,
      size_t vertex_count_limit = 10,
      size_t max_coverage = std::numeric_limits<size_t>::max())
      : base(g),
        min_coverage_gap_(min_coverage_gap),
        local_coverage_f_(local_coverage_f),
        vertex_count_limit_(vertex_count_limit),
        max_coverage_(max_coverage),
        component_remover_(g) {
  }

 protected:

  /*virtual*/
  bool ProcessEdge(EdgeId e) {

  }

};

}

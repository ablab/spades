//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <limits>

#include "visualization/visualization.hpp"
#include "compressor.hpp"
#include "dominated_set_finder.hpp"
#include "assembly_graph/graph_support/parallel_processing.hpp"


namespace omnigraph{

template<class Graph>
class ComplexTipFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;

    double relative_coverage_treshold_;
    size_t edge_length_treshold_;
    size_t max_path_length_;

    double GetTipCoverage(const GraphComponent<Graph>& component) const {
        double cov = numeric_limits<double>::max();
        for (auto edge : component.edges()) {
            cov = std::min(cov, g_.coverage(edge));
        }
        return cov;
    }

    double GetOutwardCoverage(const GraphComponent<Graph>& component) const {
        double cov = 0.0;
        for (auto v : component.vertices()) {
            for (auto edge : g_.IncidentEdges(v)) {
                if (!component.contains(edge)) {
                    cov = std::max(cov, g_.coverage(edge));
                }
            }
        }
        return cov;
    }

    double GetRelativeTipCoverage(const GraphComponent<Graph>& component) const {
        return GetTipCoverage(component) / GetOutwardCoverage(component);
    }

    bool ComponentCheck(const GraphComponent<Graph>& component) const {
        if (component.empty() || component.e_size() == 0)
            return false;

        //check if usual tip
        if (component.vertices().size() == 2) {
            DEBUG("Component is a tip! Exiting...");
            return false;
        }

        //checking edge lengths
        if (std::any_of(component.e_begin(), component.e_end(), [&](EdgeId e) {return g_.length(e) > edge_length_treshold_;})) {
            DEBUG("Tip contains too long edges");
            return false;
        }

        if (math::ge(GetRelativeTipCoverage(component), relative_coverage_treshold_)) {
            DEBUG("Tip is too high covered with respect to external edges");
            return false;
        }

        return true;
    }

public:
    ComplexTipFinder(const Graph& g, double relative_coverage,
                      size_t max_edge_length, size_t max_path_length)
            : g_(g),
              relative_coverage_treshold_(math::ge(relative_coverage, 0.0) ?
                                          relative_coverage : std::numeric_limits<double>::max()),
              edge_length_treshold_(max_edge_length), max_path_length_(max_path_length)
    { }

    GraphComponent<Graph> operator()(VertexId v) const {
        GraphComponent<Graph> empty(g_);
        VERIFY(empty.empty());
        if (g_.IncomingEdgeCount(v) != 0) {
            return empty;
        }

        DominatedSetFinder<Graph> finder(g_, v, max_path_length_);
        if (finder.FillDominated()) {
            auto ranges = finder.dominated();
            auto dom_component = finder.AsGraphComponent();
            std::set<EdgeId> component_edges(dom_component.edges());
            for (auto v : dom_component.exits()) {
                size_t current_path_length = ranges[v].end_pos;
                for (auto e : g_.OutgoingEdges(v)) {
                    if (current_path_length + g_.length(e) > max_path_length_) {
                        DEBUG("Component contains too long paths");
                        return empty;
                    }
                    component_edges.insert(e);
                }
            }
            auto extended_component = GraphComponent<Graph>::FromEdges(g_, component_edges);
            if (ComponentCheck(extended_component))
                return extended_component;
            else
                return empty;
        } else {
            DEBUG("Failed to find dominated component");
            return empty;
        }
    }

private:
    DECL_LOGGER("ComplexTipClipper")
};

template<class Graph>
class ComplexTipClipper : public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentProcessingAlgorithm<Graph, VertexId> base;
    typedef typename ComponentRemover<Graph>::HandlerF HandlerF;

    string pics_folder_;
    ComplexTipFinder<Graph> finder_;
    ComponentRemover<Graph> component_remover_;

public:
    //track_changes=false leads to every iteration run from scratch
    ComplexTipClipper(Graph& g, double relative_coverage,
                      size_t max_edge_len, size_t max_path_len,
                      size_t chunk_cnt,
                      const string& pics_folder = "" ,
                      HandlerF removal_handler = nullptr) :
            base(g, nullptr, false, std::less<VertexId>(), /*track changes*/false),
            pics_folder_(pics_folder),
            finder_(g, relative_coverage, max_edge_len, max_path_len),
            component_remover_(g, removal_handler) {
        if (!pics_folder_.empty()) {
            make_dir(pics_folder_);
        }
        this->interest_el_finder_ = std::make_shared<ParallelInterestingElementFinder<Graph, VertexId>>(
                [&](VertexId v) {return !finder_(v).empty();}, chunk_cnt);
    }

    bool Process(VertexId v) override {
        DEBUG("Processing vertex " << this->g().str(v));
        auto component = finder_(v);
        if (component.empty()) {
            DEBUG("Failed to detect complex tip starting with vertex " << this->g().str(v));
            return false;
        }

        if (!pics_folder_.empty()) {
            visualization::visualization_utils::WriteComponentSinksSources(component,
                                                      pics_folder_
                                                      + std::to_string(this->g().int_id(v)) //+ "_" + std::to_string(candidate_cnt)
                                                      + ".dot");
        }

        VERIFY(component.e_size() && component.v_size());
        DEBUG("Detected tip component edge cnt: " << component.e_size());
        component_remover_.DeleteComponent(component.e_begin(), component.e_end());
        DEBUG("Complex tip removed");
        return true;
    }

private:
    DECL_LOGGER("ComplexTipClipper")
};

}

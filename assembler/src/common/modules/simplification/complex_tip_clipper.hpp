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


namespace omnigraph{


template<class Graph>
class ComplexTipSetFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    size_t max_length_;
    std::map<VertexId, Range> complex_tc_vertices_;
    DECL_LOGGER("ComplexTipClipper")
public:
    ComplexTipSetFinder(const Graph& g, size_t max_length = std::numeric_limits<size_t>::max())
            : g_(g), max_length_(max_length)
    {    }


    GraphComponent<Graph> GetComponent(VertexId v) {
        if(g_.IncomingEdgeCount(v) != 0) {
            return GraphComponent<Graph>(g_, false);
        }
        DominatedSetFinder<Graph> finder(g_, v, max_length_);
        if(finder.FillDominated()) {
            auto ranges = finder.dominated();
            GraphComponent<Graph> dom_component = finder.AsGraphComponent();
            GraphComponent<Graph> result(dom_component);
            set<EdgeId> edges_to_add;
            for(auto v : dom_component.sinks()) {
                size_t current_path_length = ranges[v].end_pos;
                for(auto e : g_.OutgoingEdges(v)) {
                    if(current_path_length + g_.length(e) > max_length_) {
                        DEBUG("Component contains too long paths");
                        return GraphComponent<Graph>(g_, false);
                    }
                    result.AddEdge(e);
                }
            }
            return result;
        } else {
            DEBUG("Component contains too long paths");
            return GraphComponent<Graph>(g_, false);
        }
    }


};


template<class Graph>
class ComplexTipClipper {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    Graph& g_;
    double relative_coverage_treshold_;
    size_t edge_length_treshold_;
    size_t max_path_length_;
    string pics_folder_;
    std::function<void(const set<EdgeId>&)> removal_handler_;

    bool CheckEdgeLenghts(const GraphComponent<Graph>& component) const {
        for(auto e : component.edges()) {
            if(g_.length(e) > edge_length_treshold_) {
                return false;
            }
        }
        return true;
    }


    bool IsCommonTip(const GraphComponent<Graph> & component) const {
        return component.vertices().size() == 2;

    }

    template<class ElemType>
    void InsertIfNotConjugate(std::set<ElemType>& elems, ElemType elem) const {
        if (elems.count(g_.conjugate(elem)) == 0) {
            elems.insert(elem);
        }
    }


    bool RemoveComplexTip(const GraphComponent<Graph>& component) {
        ComponentRemover<Graph> remover(g_, removal_handler_);
        remover.DeleteComponent(component.edges().begin(), component.edges().end());
        return true;
    }


    bool CheckPathLengths(const map<VertexId, Range>& ranges) const {
        for(auto r : ranges) {
            if(r.second.start_pos > max_path_length_) {
                return false;
            }
        }
        return true;
    }

    double GetTipCoverage(const GraphComponent<Graph> & component) const {
        double cov = numeric_limits<double>::max();
        for(auto edge : component.edges()) {
            cov = std::min(cov, g_.coverage(edge));
        }
        return cov;
    }

    double GetOutwardCoverage(const GraphComponent<Graph> & component) const {
        double cov = 0.0;
        for(auto v : component.vertices()) {
            for(auto edge : g_.OutgoingEdges(v)) {
                if(component.contains(edge)) {
                    cov = max(cov, g_.coverage(edge));
                }
            }

            for(auto edge : g_.IncomingEdges(v)) {
                if(component.contains(edge)) {
                    cov = max(cov, g_.coverage(edge));
                }
            }
        }
        return cov;
    }

    double GetRelativeTipCoverage(const GraphComponent<Graph> & component) const {
        return GetTipCoverage(component) / GetOutwardCoverage(component);
    }

public:
    ComplexTipClipper(Graph& g, double relative_coverage, size_t max_edge_len, size_t max_path_len, const string& pics_folder = "" , std::function<void(const set<EdgeId>&)> removal_handler = 0) :
            g_(g), relative_coverage_treshold_(math::ge(relative_coverage, 0.0) ? relative_coverage : std::numeric_limits<double>::max()), edge_length_treshold_(max_edge_len) ,max_path_length_(max_path_len), pics_folder_(pics_folder), removal_handler_(removal_handler)
    { }

    bool Run() {
        size_t cnt = 0;
        INFO("Complex tip clipper started");
        if (!pics_folder_.empty()) {
            make_dir(pics_folder_);
        }

        bool something_done_flag = false;
        ComplexTipSetFinder<Graph> complex_tc_set_finder(g_, max_path_length_);
        for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
            if(g_.IncomingEdgeCount(*it) != 0) {
                continue;
            }
            DEBUG("Processing vertex " << g_.str(*it));

            GraphComponent<Graph> component = complex_tc_set_finder.GetComponent(*it);
            if(component.v_size() == 0) {
                continue;
            }

            if(!CheckEdgeLenghts(component)) {
                DEBUG("Tip contains too long edges");
                continue;
            }

            if(IsCommonTip(component)) {
                DEBUG("Component is a tip! Exiting...");
                continue;
            }


            if(math::ge(GetRelativeTipCoverage(component), relative_coverage_treshold_)) {
                DEBUG("Tip is too high covered with respect to external edges");
                continue;
            }

            if (!pics_folder_.empty()) {
                visualization::WriteComponentSinksSources(component,
                        pics_folder_
                                + ToString(g_.int_id(*it)) //+ "_" + ToString(candidate_cnt)
                                + ".dot");
            }
            DEBUG(component.vertices().size());
            something_done_flag = RemoveComplexTip(component);

            if(something_done_flag) {
                cnt++;
                DEBUG("Removed");
            }
        }

        CompressAllVertices(g_);
        DEBUG("Complex tip clipper finished");
        DEBUG("Tips processed " << cnt);
        return something_done_flag;
    }
private:
    DECL_LOGGER("ComplexTipClipper")
};

}

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


    bool CheckSize(const GraphComponent<Graph> & component) const {
        return (component.vertices().size() > 1);
    }

    void RemoveComplexTip(GraphComponent<Graph>& component) {
        ComponentRemover<Graph> remover(g_, removal_handler_);
        remover.DeleteComponent(component.edges().begin(), component.edges().end());
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
    ComplexTipClipper(Graph& g, double relative_coverage, size_t max_edge_len, size_t max_path_len, const string& pics_folder = "", std::function<void(const set<EdgeId>&)> removal_handler = 0) :
            g_(g), relative_coverage_treshold_(math::ge(relative_coverage, 0.0) ? relative_coverage : std::numeric_limits<double>::max()), edge_length_treshold_(max_edge_len) ,max_path_length_(max_path_len), pics_folder_(pics_folder), removal_handler_(removal_handler)
    { }

    bool Run() {
        size_t cnt = 0;
        INFO("Complex tip clipper started");
        if (!pics_folder_.empty()) {
            make_dir(pics_folder_);
        }

        bool something_done_flag = false;
        for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
            if(g_.IncomingEdgeCount(*it) != 0) {
                continue;
            }
            DEBUG("Processing vertex " << g_.str(*it));

            DominatedSetFinder<Graph> dom_finder(g_, *it, max_path_length_ * 2);
            dom_finder.FillDominated();
            auto component = dom_finder.AsGraphComponent();

            if(!CheckEdgeLenghts(component)) {
                DEBUG("Tip contains too long edges");
                continue;
            }

            if(!CheckSize(component)) {
                DEBUG("Component doesn't meet size requirements");
                continue;
            }
            auto dominated = dom_finder.dominated();
            if(!CheckPathLengths(dominated)) {
                DEBUG("Tip contains too long paths");
                continue;
            }

            if(math::ls(GetRelativeTipCoverage(component), relative_coverage_treshold_)) {
                DEBUG("Tip is too high covered with respect to external edges");
                continue;
            }

            if (!pics_folder_.empty()) {
                visualization::WriteComponentSinksSources(component,
                        pics_folder_
                                + ToString(g_.int_id(*it)) //+ "_" + ToString(candidate_cnt)
                                + ".dot");
            }

            something_done_flag = true;
            cnt++;
            RemoveComplexTip(component);
        }
        CompressAllVertices(g_);
        INFO("Complex tip clipper finished");
        INFO("Tips processed " << cnt);
        return something_done_flag;
    }
private:
    DECL_LOGGER("ComplexTipClipper")
};

}

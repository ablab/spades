#pragma once

#include <limits>

#include "omni_utils.hpp"
#include "omni/visualization/visualization.hpp"


namespace omnigraph{


template<class Graph>
class ComplexTipClipper {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    Graph& g_;
    size_t max_length_;
    string pics_folder_;

    const size_t edge_length_treshold = 100;

    bool CheckEdgeLenghts(const GraphComponent<Graph>& component) const {
        for(auto e : component.edges()) {
            if(g_.length(e) > edge_length_treshold) {
                return false;
            }
        }
        return true;
    }


    bool CheckSize(const GraphComponent<Graph> & component) const {
        return (component.vertices().size() > 1);
    }

    void RemoveComplexTip(GraphComponent<Graph>& component) {
        ComponentRemover<Graph> remover(g_);
        remover.DeleteComponent(component.edges().begin(), component.edges().end());
    }


    bool CheckPathLengths(const map<VertexId, Range>& ranges) const {
        for(auto r : ranges) {
            if(r.second.start_pos > max_length_) {
                return false;
            }
        }
        return true;
    }

public:
    ComplexTipClipper(Graph& g, size_t max_length, const string& pics_folder = "") :
            g_(g), max_length_(max_length), pics_folder_(pics_folder)
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

            {
                DominatedSetFinder<Graph> dom_finder(g_, *it, max_length_ * 2);
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

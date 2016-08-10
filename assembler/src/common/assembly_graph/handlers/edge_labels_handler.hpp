//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 *
 * Saves labeling of new_graph via different graph transformation by edges of unresolved graph - old_graph
 * Has two methods
 *
 *  Created on: Aug 5, 2011
 *      Author: undead
 */

#ifndef EDGE_LABELS_HANDLER_HPP_
#define EDGE_LABELS_HANDLER_HPP_

//#include "utils.hpp"
#include "visualization/graph_labeler.hpp"
#include "utils/simple_tools.hpp"
#include <unordered_map>
#include <map>

using namespace omnigraph;

namespace omnigraph {
using std::map;

//todo ask Shurik to remove new_graph_
template<class Graph>
class EdgeLabelHandler : public GraphActionHandler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
private:
    Graph &new_graph_;
    Graph &old_graph_;
    //From new edge to sequence of old
public:
    map<EdgeId, vector<EdgeId> > edge_labels;
    //From old edge to set of new ones, containing it.
    map<EdgeId, set<EdgeId> > edge_inclusions;
public:
    //TODO: integrate this to resolver, remove "from_resolve" parameter
    EdgeLabelHandler(Graph &new_graph, Graph &old_graph,
                     const std::map<EdgeId, EdgeId> &from_resolve)
            : GraphActionHandler<Graph>(new_graph, "EdgePositionHandler"),
              new_graph_(new_graph),
              old_graph_(old_graph) {
        // printing from resolve
        FillLabels(from_resolve);
        /*        for(auto iter = from_resolve.begin(); iter != from_resolve.end(); ++iter) {
         if (edge_inclusions.find(iter->second) == edge_inclusions.end()){
         set<EdgeId> tmp;
         edge_inclusions.insert(make_pair(iter->second, tmp));
         }
         edge_inclusions[iter->second].insert(iter->first);

         if (edge_labels.find(iter->first) == edge_labels.end()) {
         set<EdgeId> tmp;
         edge_labels.insert(make_pair(iter->first, tmp));
         }
         edge_labels[iter->second].push_back(iter->second);
         }
         */}

    EdgeLabelHandler(Graph &new_graph, Graph &old_graph)
            : GraphActionHandler<Graph>(new_graph, "EdgePositionHandler"),
              new_graph_(new_graph),
              old_graph_(old_graph) {
    }

    void FillLabels(const map<EdgeId, EdgeId> &from_resolve) {
        for (auto iter = from_resolve.begin(); iter != from_resolve.end();
             ++iter) {
            if (edge_inclusions.find(iter->second) == edge_inclusions.end()) {
                set<EdgeId> tmp;
                edge_inclusions.insert(make_pair(iter->second, tmp));
            }
            edge_inclusions.find(iter->second)->second.insert(iter->first);

            if (edge_labels.find(iter->first) == edge_labels.end()) {
                vector<EdgeId> tmp;
                edge_labels.insert(make_pair(iter->first, tmp));
            }
            edge_labels[iter->first].push_back(iter->second);
        }
    }

    virtual ~EdgeLabelHandler() {
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
        TRACE("Handle glue");
        if (edge_labels[edge1] != edge_labels[edge2])
            WARN("gluing two different edges is not a good idea on this step! EdgeLabel Handler can fail on such operation");
        vector<EdgeId> tmp;
        for (size_t i = 0; i < edge_labels[edge1].size(); i++) {
            edge_inclusions.find(edge_labels[edge1][i])->second.insert(
                    new_edge);
            edge_inclusions.find(edge_labels[edge1][i])->second.erase(edge1);
            tmp.push_back(edge_labels[edge1][i]);

            edge_labels.erase(edge1);
        }
        for (size_t i = 0; i < edge_labels[edge2].size(); i++) {
            edge_inclusions.find(edge_labels[edge2][i])->second.insert(
                    new_edge);
            edge_inclusions.find(edge_labels[edge2][i])->second.erase(edge2);
            edge_labels.erase(edge2);

            //    tmp.push_back(edge_labels[edge1][i]);
        }

        edge_labels.insert(make_pair(new_edge, tmp));

    }

    virtual void HandleSplit(EdgeId /*oldEdge*/, EdgeId /*newEdge1*/, EdgeId /*newEdge2*/) {
        WARN("EdgesLabelHandler does not support splits");
    }

    virtual void HandleMerge(const vector<EdgeId> &oldEdges, EdgeId newEdge) {
        TRACE("HandleMerge by edge labels handler");
        size_t n = oldEdges.size();
        vector<EdgeId> tmp;
        for (size_t j = 0; j < n; j++) {
            TRACE("Edge " << oldEdges[j] << " was labeled by " << edge_labels[oldEdges[j]]);
            for (size_t i = 0; i < edge_labels[oldEdges[j]].size(); i++) {
                edge_inclusions[edge_labels[oldEdges[j]][i]].insert(newEdge);
                edge_inclusions[edge_labels[oldEdges[j]][i]].erase(oldEdges[j]);
                tmp.push_back(edge_labels[oldEdges[j]][i]);
            }
            edge_labels.erase(oldEdges[j]);
        }
        if (edge_labels.find(newEdge) != edge_labels.end()) {
            DEBUG("Unexpected finding of new edge labels");
        };
        edge_labels[newEdge] = tmp;

    }

    /*
     virtual void HandleAdd(VertexId v) {
     AddVertexIntId(v);
     }
     virtual void HandleDelete(VertexId v) {
     ClearVertexId(v);
     }
     */
    virtual void HandleAdd(EdgeId e) {
        TRACE("Add edge " << e);

    }

    virtual void HandleDelete(EdgeId e) {
        for (size_t i = 0; i < edge_labels[e].size(); i++) {
            edge_inclusions[edge_labels[e][i]].erase(e);
        }
        edge_labels.erase(e);
    }

    std::string str(EdgeId edgeId) const {
        std::stringstream ss;

        auto it = edge_labels.find(edgeId);
        if (it != edge_labels.end()) {
            TRACE("Number of labels " << it->second.size());
            for (auto label_it = it->second.begin(), end = it->second.end();
                 label_it != end; ++label_it) {
                ss << this->g().str(*label_it) << "\\n";
            }
        }
        return ss.str();
    }

    vector<pair<EdgeId, size_t> > resolvedPositions(EdgeId old_edge, size_t position_on_edge) {
        vector<pair<EdgeId, size_t> > res;
        for (auto it = edge_inclusions[old_edge].begin(); it != edge_inclusions[old_edge].end(); it++) {
            EdgeId cur_edge = *it;
            size_t cur_shift = 0;
            for (size_t i = 0; i < edge_labels[cur_edge].size(); i++) {
                if (edge_labels[cur_edge][i] == old_edge) {
                    res.push_back(make_pair(cur_edge, cur_shift + position_on_edge));
                }
                cur_shift += old_graph_.length(edge_labels[cur_edge][i]);
            }
        }
        return res;
    }

};

template<class Graph>
class EdgesLabelsGraphLabeler : public GraphLabeler<Graph> {

protected:
    typedef GraphLabeler<Graph> super;
    typedef typename super::EdgeId EdgeId;
    typedef typename super::VertexId VertexId;
    Graph &g_;
public:
    EdgeLabelHandler<Graph> &EdgesLabels;

    EdgesLabelsGraphLabeler(Graph &g, EdgeLabelHandler<Graph> &EdgesLab)
            : g_(g),
              EdgesLabels(EdgesLab) {
    }

    virtual std::string label(VertexId vertexId) const {
        return g_.str(vertexId);
    }

    virtual std::string label(EdgeId edgeId) const {
        return EdgesLabels.str(edgeId) + ": " + g_.str(edgeId);
    }

    virtual ~EdgesLabelsGraphLabeler() {
        TRACE("~EdgesPosGraphLabeler");
    }

};
}

#endif /* EDGE_LABELS_HANDLER_HPP_ */

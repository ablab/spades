/*
 * path_visualizer.hpp
 *
 *  Created on: Mar 22, 2012
 *      Author: andrey
 */

#ifndef PATH_VISUALIZER_HPP_
#define PATH_VISUALIZER_HPP_

#include "bidirectional_path.hpp"

namespace path_extend {

using namespace debruijn_graph;

template<class Graph>
class PathGraphLabeler : public AbstractGraphLabeler<Graph> {
    typedef AbstractGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    std::map<EdgeId, std::string> labels_;

public:
    PathGraphLabeler(Graph& g, PathContainer& paths) : base(g) {
//      for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//          labels_[*iter] = "";
//      }

        for(size_t i = 0; i < paths.size(); ++i) {
            BidirectionalPath * path = paths.Get(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (labels_.count(path->At(j)) > 0) {
                    labels_[path->At(j)] += ", ";
                }
                labels_[path->At(j)] += "(" + ToString(path->GetId()) + " : " + ToString(j) + ")";
            }

            path = paths.GetConjugate(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (labels_.count(path->At(j)) > 0) {
                    labels_[path->At(j)] += ", ";
                }
                labels_[path->At(j)] += "(" + ToString(path->GetId()) + " : " + ToString(j) + ")";
            }
        }
    }

    virtual std::string label(VertexId vertexId) const {
        return "";
    }

    virtual std::string label(EdgeId edgeId) const {
        auto label = labels_.find(edgeId);
        return label == labels_.end() ? "" : label->second;
    }
};


class PathVisualizer {

protected:
    bool writeLength;
    bool writePos;

    size_t k_;

public:

    PathVisualizer(size_t k): writeLength(true), writePos(true), k_(k) {

    }

    void writeGraphWithPathsSimple(conj_graph_pack& gp, const string& file_name, const string& graph_name, PathContainer& paths) {
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        StrGraphLabeler<Graph> str_labeler(gp.g);
        PathGraphLabeler<Graph> path_labeler(gp.g, paths);
        EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

        CompositeLabeler<Graph> composite_labeler(str_labeler, path_labeler, pos_labeler);

        auto_ptr<GraphColorer<Graph>> colorer(DefaultColorer(gp.g, FindGenomePath(gp.genome, gp.g, gp.index, k_)
                , FindGenomePath(!gp.genome, gp.g, gp.index, k_)));

        omnigraph::DotGraphPrinter<Graph> printer(gp.g, composite_labeler, *colorer, graph_name, filestr);
        ColoredGraphVisualizer<Graph> gv(gp.g, printer);
        AdapterGraphVisualizer<Graph> result_vis(gp.g, gv);
        result_vis.Visualize();
        filestr.close();
        INFO("Visualizing graph " << graph_name << " done");
    }

    bool isWriteLength() const
    {
        return writeLength;
    }

    bool isWritePos() const
    {
        return writePos;
    }

    void setWriteLength(bool writeLength)
    {
        this->writeLength = writeLength;
    }

    void setWritePos(bool writePos)
    {
        this->writePos = writePos;
    }
};

}

#endif /* PATH_VISUALIZER_HPP_ */

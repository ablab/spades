//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_visualizer.hpp
 *
 *  Created on: Mar 22, 2012
 *      Author: andrey
 */

#ifndef PATH_VISUALIZER_HPP_
#define PATH_VISUALIZER_HPP_

#include "bidirectional_path.hpp"
#include "stats/debruijn_stats.hpp"

namespace path_extend {

using namespace debruijn_graph;

template<class Graph>
class PathGraphLabeler : public AbstractGraphLabeler<Graph> {
    typedef AbstractGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    std::map<EdgeId, std::string> labels_;

public:
    PathGraphLabeler(const Graph& g, PathContainer& paths) : base(g) {
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

    virtual std::string label(VertexId /*vertexId*/) const {
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

public:

    PathVisualizer(): writeLength(true), writePos(true) {

    }

    void writeGraphWithPathsSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name, PathContainer& paths) const {
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        StrGraphLabeler<Graph> str_labeler(gp.g);
        PathGraphLabeler<Graph> path_labeler(gp.g, paths);
        CoverageGraphLabeler<Graph> cov_labler(gp.g);
        EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

        CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, path_labeler, pos_labeler);
        shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer;
        if (gp.index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            colorer = omnigraph::visualization::DefaultColorer(gp.g);
        }

        omnigraph::visualization::ComponentVisualizer<Graph> visualizer(gp.g, false);
        omnigraph::visualization::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        StrGraphLabeler<Graph> str_labeler(gp.g);
        EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
        CoverageGraphLabeler<Graph> cov_labler(gp.g);
        CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, pos_labeler);

        shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer;

        if (gp.index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            Path<EdgeId> empty;
            colorer = omnigraph::visualization::DefaultColorer(gp.g, empty, empty);
        }

        omnigraph::visualization::ComponentVisualizer<Graph> visualizer(gp.g, false);
        omnigraph::visualization::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const Graph& g, const string& file_name, const string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        StrGraphLabeler<Graph> str_labeler(g);
        CoverageGraphLabeler<Graph> cov_labler(g);
        CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler);

        shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer;

        Path<EdgeId> empty;
        colorer = omnigraph::visualization::DefaultColorer(g, empty, empty);

        omnigraph::visualization::ComponentVisualizer<Graph> visualizer(g, false);
        omnigraph::visualization::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
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

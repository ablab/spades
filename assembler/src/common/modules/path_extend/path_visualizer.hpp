//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * path_visualizer.hpp
 *
 *  Created on: Mar 22, 2012
 *      Author: andrey
 */

#ifndef PATH_VISUALIZER_HPP_
#define PATH_VISUALIZER_HPP_

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

namespace path_extend {

using namespace debruijn_graph;

template<class Graph>
class PathGraphLabeler : public visualization::graph_labeler::AbstractGraphLabeler<Graph> {
    typedef visualization::graph_labeler::AbstractGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    std::map<EdgeId, std::string> labels_;

public:
    PathGraphLabeler(const Graph& g, const PathContainer& paths) : base(g) {
        for(size_t i = 0; i < paths.size(); ++i) {
            BidirectionalPath * path = paths.Get(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (labels_.count(path->At(j)) > 0) {
                    labels_[path->At(j)] += ", ";
                }
                labels_[path->At(j)] += "(" + std::to_string(path->GetId()) + " : " + std::to_string(j) + ")";
            }

            path = paths.GetConjugate(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (labels_.count(path->At(j)) > 0) {
                    labels_[path->At(j)] += ", ";
                }
                labels_[path->At(j)] += "(" + std::to_string(path->GetId()) + " : " + std::to_string(j) + ")";
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

    void writeGraphWithPathsSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name, const PathContainer& paths) const {
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(gp.g);
        PathGraphLabeler<Graph> path_labeler(gp.g, paths);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(gp.g);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);

        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, path_labeler, pos_labeler);
        shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;
        if (gp.index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            colorer = visualization::graph_colorer::DefaultColorer(gp.g);
        }

        visualization::visualizers::ComponentVisualizer<Graph> visualizer(gp.g, false);
        visualization::vertex_linker::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const conj_graph_pack& gp, const string& file_name, const string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(gp.g);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(gp.g);
        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, pos_labeler);

        shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;

        if (gp.index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            Path<EdgeId> empty;
            colorer = visualization::graph_colorer::DefaultColorer(gp.g, empty, empty);
        }

        visualization::visualizers::ComponentVisualizer<Graph> visualizer(gp.g, false);
        visualization::vertex_linker::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const Graph& g, const string& file_name, const string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(g);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(g);
        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler);

        shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;

        Path<EdgeId> empty;
        colorer = visualization::graph_colorer::DefaultColorer(g, empty, empty);

        visualization::visualizers::ComponentVisualizer<Graph> visualizer(g, false);
        visualization::vertex_linker::EmptyGraphLinker<Graph> linker;
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

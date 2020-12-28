//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef PATH_VISUALIZER_HPP_
#define PATH_VISUALIZER_HPP_

#include "assembly_graph/paths/bidirectional_path_container.hpp"
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
            {
                const BidirectionalPath &path = paths.Get(i);
                for (size_t j = 0; j < path.Size(); ++j) {
                    if (labels_.count(path.At(j)) > 0) {
                        labels_[path.At(j)] += ", ";
                    }
                    labels_[path.At(j)] += "(" + std::to_string(path.GetId()) + " : " + std::to_string(j) + ")";
                }
            }
            

            {
                const BidirectionalPath &path = paths.GetConjugate(i);
                for (size_t j = 0; j < path.Size(); ++j) {
                    if (labels_.count(path.At(j)) > 0) {
                        labels_[path.At(j)] += ", ";
                    }
                    labels_[path.At(j)] += "(" + std::to_string(path.GetId()) + " : " + std::to_string(j) + ")";
                }
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


//TODO: refactor this copy-pasta
class PathVisualizer {

protected:
    bool writeLength;
    bool writePos;

public:

    PathVisualizer(): writeLength(true), writePos(true) {

    }

    void writeGraphWithPathsSimple(const GraphPack& gp, const std::string& file_name,
                                   const std::string& graph_name, const PathContainer& paths) const {
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        const auto &graph = gp.get<Graph>();
        const auto &edge_pos = gp.get<EdgesPositionHandler<Graph>>();
        const auto &index = gp.get<EdgeIndex<Graph>>();

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(graph);
        PathGraphLabeler<Graph> path_labeler(graph, paths);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(graph);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(graph, edge_pos);

        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, path_labeler, pos_labeler);
        std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;
        if (index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            colorer = visualization::graph_colorer::DefaultColorer(graph);
        }

        visualization::visualizers::ComponentVisualizer<Graph> visualizer(graph, false);
        visualization::vertex_linker::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const GraphPack& gp, const std::string& file_name, const std::string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        const auto &graph = gp.get<Graph>();
        const auto &edge_pos = gp.get<EdgesPositionHandler<Graph>>();
        const auto &index = gp.get<EdgeIndex<Graph>>();

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(graph);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(graph, edge_pos);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(graph);
        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler, pos_labeler);

        std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;

        if (index.IsAttached()) {
             colorer = stats::DefaultColorer(gp);
        } else {
            Path<EdgeId> empty;
            colorer = visualization::graph_colorer::DefaultColorer(graph, empty, empty);
        }

        visualization::visualizers::ComponentVisualizer<Graph> visualizer(graph, false);
        visualization::vertex_linker::EmptyGraphLinker<Graph> linker;
        visualizer.Visualize(filestr, composite_labeler, *colorer, linker);
        filestr.close();
        INFO("Visualizing graph done");
    }

    void writeGraphSimple(const Graph& g, const std::string& file_name, const std::string& graph_name) const{
        INFO("Visualizing graph " << graph_name << " to file " << file_name);
        std::fstream filestr;
        filestr.open(file_name.c_str(), std::fstream::out);

        visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(g);
        visualization::graph_labeler::CoverageGraphLabeler<Graph> cov_labler(g);
        visualization::graph_labeler::CompositeLabeler<Graph> composite_labeler(str_labeler, cov_labler);

        std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer;

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

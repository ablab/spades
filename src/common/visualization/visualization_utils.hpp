#pragma once

//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_printer.hpp"
#include "visualizers.hpp"
#include "vertex_linker.hpp"

#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/components/splitters.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "visualizers.hpp"
#include "vertex_linker.hpp"

#include <filesystem>
#include <fstream>

namespace visualization {

namespace visualization_utils {

template<class Graph>
void WriteComponents(const Graph &g,
                     const std::filesystem::path &folder_name,
                     std::shared_ptr<GraphSplitter<Graph>> inner_splitter,
                     std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                     const graph_labeler::GraphLabeler<Graph> &labeler) {
    vertex_linker::EmptyGraphLinker<Graph> linker;
//  shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
    auto filter = std::make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    std::shared_ptr<GraphSplitter<Graph>> splitter = std::make_shared<omnigraph::CollectingSplitterWrapper<Graph>>(
            inner_splitter, filter);
    visualization::visualizers::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(*splitter,
                                                                                                   folder_name);
}

template<class Graph>
void DrawComponentsOfShortEdges(const Graph &g, const std::filesystem::path &output_dir, size_t min_length, size_t sinks,
                                size_t sources) {
    std::vector<typename Graph::EdgeId> short_edges;
    std::filesystem::path pics_folder_ =
            output_dir / (std::to_string(min_length) + "_" + std::to_string(sinks) + "_" + std::to_string(sources) + "_" + "pics_polymorphic/");
    create_directory(pics_folder_);
    INFO("Writing pics with components consisting of short edges to " + static_cast<std::string>(pics_folder_));
    auto splitter = LongEdgesExclusiveSplitter<Graph>(g, min_length);
    while (splitter->HasNext()) {
        GraphComponent<Graph> component = splitter->Next();
        if (component.v_size() > 3 && component.exits().size() == sinks &&
                component.entrances().size() == sources) {
            bool fail = false;
            for (auto v : component.entrances()) {
                if (component.g().IncomingEdgeCount(v) != 1) {
                    fail = true;
                }
            }
            for (auto v : component.exits()) {
                if (component.g().OutgoingEdgeCount(v) != 1) {
                    fail = true;
                }
            }

            if (fail) {
                continue;
            }

            graph_labeler::StrGraphLabeler<Graph> labeler(component.g());
            graph_labeler::CoverageGraphLabeler<Graph> labeler2(component.g());
            graph_labeler::CompositeLabeler<Graph> compositeLabeler(labeler, labeler2);
            WriteComponentSinksSources(component,
                                       pics_folder_ + std::to_string(g.int_id(*component.vertices().begin()))
                                       + ".dot", visualization::graph_colorer::DefaultColorer(g),
                                       compositeLabeler);
            INFO("Component is written to " + std::to_string(g.int_id(*component.vertices().begin())) + ".dot");

            //            PrintComponent(component,
//                                pics_folder_ + "ShortComponents/"
//                                        + std::to_string(gp.g.int_id(component.vertices_[0]))
//                                         + ".dot");
        }
    }
}


template<class Graph>
void WriteSizeLimitedComponents(const Graph &g,
                                const std::filesystem::path &folder_name,
                                std::shared_ptr<GraphSplitter<Graph>> inner_splitter,
                                std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                                const graph_labeler::GraphLabeler<Graph> &labeler, int min_component_size,
                                int max_component_size, size_t max_components) {
    vertex_linker::EmptyGraphLinker<Graph> linker;

    auto filter = std::make_shared<omnigraph::ComponentSizeFilter<Graph>>(g, 1000000000,
                                                                          (size_t) min_component_size,
                                                                          (size_t) max_component_size);
    auto splitter = std::make_shared<omnigraph::CollectingSplitterWrapper<Graph>>(inner_splitter, filter);
    visualization::visualizers::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker, false,
                                                                max_components).SplitAndVisualize(*splitter, folder_name);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph> &gc,
                    const std::filesystem::path &file_name, std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                    const graph_labeler::GraphLabeler<Graph> &labeler) {
    vertex_linker::EmptyGraphLinker<Graph> linker;
    graph_colorer::BorderDecorator<Graph> component_colorer(gc, *colorer, "yellow");
    std::ofstream os;
    os.open(file_name);
    visualization::visualizers::ComponentVisualizer<Graph>(gc.g(), true).
            Visualize(gc, os, labeler, component_colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentSinksSources(const GraphComponent<Graph> &gc,
                                const std::filesystem::path &file_name, std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                                const graph_labeler::GraphLabeler<Graph> &labeler) {
    vertex_linker::EmptyGraphLinker<Graph> linker;
    graph_colorer::SinkSourceDecorator<Graph> component_colorer(gc, *colorer);
    std::ofstream os;
    os.open(file_name);
    visualization::visualizers::ComponentVisualizer<Graph>(gc.g(), true).
            Visualize(gc, os, labeler, component_colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentSinksSources(const GraphComponent<Graph> &gc,
                                const std::filesystem::path &file_name) {

    graph_labeler::StrGraphLabeler<Graph> labeler(gc.g());
    graph_labeler::CoverageGraphLabeler<Graph> labeler2(gc.g());
    graph_labeler::CompositeLabeler<Graph> compositeLabeler(labeler, labeler2);
    vertex_linker::EmptyGraphLinker<Graph> linker;
    WriteComponentSinksSources(gc, file_name, graph_colorer::DefaultColorer(gc.g()),
                               compositeLabeler);
}

template<class Graph>
void WriteSimpleComponent(const GraphComponent<Graph> &gc, const std::filesystem::path &file_name,
                          std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                          const graph_labeler::GraphLabeler<Graph> &labeler) {
    vertex_linker::EmptyGraphLinker<Graph> linker;
    std::ofstream os;
    os.open(file_name);
    visualization::visualizers::ComponentVisualizer<Graph>(gc.g(), false).
            Visualize(gc, os, labeler, *colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentsAlongPath(const Graph &g, const std::vector<typename Graph::EdgeId> &path,
                              const std::filesystem::path &prefix_path, std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer,
                              const graph_labeler::GraphLabeler<Graph> &labeler, bool color_path = true) {
    auto edge_colorer = std::make_shared<graph_colorer::CompositeEdgeColorer<Graph>>("black");
    edge_colorer->AddColorer(colorer);
    if (color_path) {
        edge_colorer->AddColorer(std::make_shared<graph_colorer::SetColorer<Graph>>(g, path, "green"));
    }
    auto resulting_colorer = std::make_shared<graph_colorer::CompositeGraphColorer<Graph>>(colorer, edge_colorer);
    auto rs = ReliableSplitterAlongPath<Graph>(g, path);
    auto filter = std::make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    auto splitter = std::make_shared<omnigraph::CondensingSplitterWrapper<Graph>>(rs, filter);
    WriteComponents<Graph>(g, prefix_path, splitter, resulting_colorer, labeler);
}

template<class Graph>
class LocalityPrintingRH {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &g_;
    const graph_labeler::GraphLabeler<Graph> &labeler_;
    std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer_;
    const std::filesystem::path output_folder_;
public:
    LocalityPrintingRH(const Graph &g, const graph_labeler::GraphLabeler<Graph> &labeler,
                       std::shared_ptr<graph_colorer::GraphColorer<Graph>> colorer, const std::filesystem::path &output_folder)
            :
            g_(g),
            labeler_(labeler),
            colorer_(colorer),
            output_folder_(output_folder) {
//        path::make_dirs(output_folder_);
    }

    void HandleDelete(EdgeId e, const std::string &add_label = "") {
        //todo magic constant
//          map<EdgeId, string> empty_coloring;
        auto edge_colorer = std::make_shared<graph_colorer::CompositeEdgeColorer<Graph>>("black");
        edge_colorer->AddColorer(colorer_);
        edge_colorer->AddColorer(
                std::make_shared<graph_colorer::SetColorer<Graph>>(g_, std::vector<EdgeId>(1, e), "green"));
        std::shared_ptr<graph_colorer::GraphColorer<Graph>> resulting_colorer =
                std::make_shared<graph_colorer::CompositeGraphColorer<Graph>>(colorer_, edge_colorer);

        std::filesystem::path fn = output_folder_ / ("edge_" + std::to_string(g_.int_id(e)) + add_label + ".dot");
        visualization::visualization_utils::WriteComponent(omnigraph::EdgeNeighborhood<Graph>(g_, e, 50, 250),
                                                           fn, resulting_colorer, labeler_);
    }

private:
    DECL_LOGGER("LocalityPrintingRH");
};

}

}

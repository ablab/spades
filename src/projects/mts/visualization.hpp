
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "visualization/graph_colorer.hpp"
#include "visualization/visualization_utils.hpp"

namespace debruijn_graph {

template <class Graph>
class AnnotatedGraphColorer
    : public visualization::graph_colorer::GraphColorer<Graph> {

    EdgeAnnotation annotation_;
    std::map<bin_id, std::string> color_map_;

public:
    AnnotatedGraphColorer(const EdgeAnnotation& annotation)
        : annotation_(annotation) {
        std::vector<std::string> preset_colors({"red", "blue", "yellow", "orange", "purple", "pink"});
        VERIFY(annotation_.interesting_bins().size() <= preset_colors.size());
        size_t i = 0;
        for (const auto& b_id : annotation_.interesting_bins()) {
            color_map_[b_id] = preset_colors[i];
            ++i;
        }
    }

    std::string GetValue(typename Graph::VertexId) const { return "black"; }

    std::string GetValue(typename Graph::EdgeId edge) const {
        if (annotation_.Annotation(edge).empty()) {
            return "black";
        }
        std::vector<std::string> colors;
        auto ann = annotation_.Annotation(edge);
        std::ostringstream ss;
        std::transform(ann.begin(), ann.end(), std::ostream_iterator<std::string>(ss, ":"), [&](bin_id b){
            return utils::get(color_map_, b);
        });
        return ss.str();
    }

};

void PrintColoredAnnotatedGraphAroundEdge(const graph_pack::GraphPack &gp, const EdgeId &edge,
                                          const EdgeAnnotation &annotation, const std::filesystem::path &output_filename) {
    //std::cout << output_filename << std::endl;
    const auto &graph = gp.get<Graph>();
    visualization::graph_labeler::DefaultLabeler<Graph> labeler(graph, gp.get<EdgesPositionHandler<Graph>>());
    auto colorer_ptr = std::make_shared<AnnotatedGraphColorer<Graph>>(annotation);
    GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(graph, edge, 100, 10000);
    visualization::visualization_utils::WriteComponent<Graph>(component, output_filename, colorer_ptr, labeler);
}

void PrintAnnotatedAlongPath(const graph_pack::GraphPack &gp, const std::vector<EdgeId> &path,
                             const EdgeAnnotation &annotation, const std::string &output_prefix) {
    const auto &graph = gp.get<Graph>();
    visualization::graph_labeler::DefaultLabeler<Graph> labeler(graph, gp.get<EdgesPositionHandler<Graph>>());
    auto colorer_ptr = std::make_shared<AnnotatedGraphColorer<Graph>>(annotation);
    visualization::visualization_utils::WriteComponentsAlongPath<Graph>(graph, path, output_prefix, colorer_ptr, labeler);
}

}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "modules/alignment/sequence_mapper.hpp"
#include "visualization/visualization_utils.hpp"

namespace cap {

//todo refactor all these methods

/*
 *  This filter restricts also components where all cops are borschs
 *  (all edges are colored in all colors)
 */
template<class Graph>
class ComponentSingleColorFilter: public GraphComponentFilter<Graph> {
private:
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    ColorHandler<Graph> color_handler_;
    TColorSet restricted_color_;
    size_t max_length_;
    size_t vertex_number_;

public:
    ComponentSingleColorFilter(const Graph &graph, const ColorHandler<Graph> &color_handler,
            const TColorSet &restricted_color, size_t max_length, size_t vertex_number)
        : base(graph),
          color_handler_(color_handler),
          restricted_color_(restricted_color),
          max_length_(max_length),
          vertex_number_(vertex_number) {
    }

    /*virtual*/
    bool Check(const GraphComponent<Graph> &gc) const {
        return true;
        TRACE("Check component");
        auto &component = gc.vertices();
        if (component.size() <= vertex_number_)
            return false;

        bool length_flag = false,
             color_flag = false;
//        set < VertexId > component(vertices.begin(), vertices.end());
        for (auto iterator = component.begin(); iterator != component.end();
                ++iterator) {
            for (EdgeId e : this->graph().OutgoingEdges(*iterator)) {
                if (component.count(this->graph().EdgeEnd(e)) == 1) {
                    if (this->graph().length(e) <= max_length_) {
                        length_flag = true;
                    }
                    if (color_handler_.Color(e) != restricted_color_) {
                        TRACE("Found good color " << color_handler_.Color(e).ToString());
                        color_flag = true;
                    }
                    if (length_flag && color_flag) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

template<class Graph>
void PrintColoredGraph(const Graph& g, const ColorHandler<Graph>& coloring,
        const EdgesPositionHandler<Graph>& pos, const string& output_filename) {
    shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitter<Graph>(g, 1000000, 30);
    visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g);
    visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

    visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
    WriteComponents(g, splitter, output_filename,
//                *ConstructColorer(coloring),
            *ConstructBorderColorer(g, coloring), labeler);
}

template<class Graph>
void PrintColoredGraphAroundEdge(const Graph& g,
    const ColorHandler<Graph>& coloring, const EdgeId edge,
    const EdgesPositionHandler<Graph>& pos, const string& output_filename) {
  INFO(output_filename);
    visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g);
    visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

    visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
    GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(g, edge);
    visualization::visualization_utils::WriteComponent(component, output_filename, coloring.ConstructColorer(component), labeler);
}

template<class Graph>
void PrintColoredGraphWithColorFilter(const Graph &g, const ColorHandler<Graph> &coloring,
    const CoordinatesHandler<Graph> &pos, const vector<string> &genome_names, const string &output_folder) {

  size_t edge_length_bound = 1000000;
  size_t colors_number = coloring.max_colors();
  TColorSet restricted_color = TColorSet::AllColorsSet(colors_number);

    shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitter<Graph>(g, edge_length_bound, 30);
    shared_ptr<omnigraph::GraphComponentFilter<Graph>> filter = make_shared<ComponentSingleColorFilter<Graph>>(g, coloring, restricted_color, edge_length_bound, 2);
    shared_ptr<omnigraph::GraphSplitter<Graph>> fs = make_shared<omnigraph::FilteringSplitterWrapper<Graph> >(splitter, filter);
    LengthIdGraphLabeler<Graph> basic_labeler(g);
    EdgeCoordinatesGraphLabeler<Graph> pos_labeler(g, pos, genome_names);

    visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
    visualization::visualization_utils::WriteComponents(g, output_folder, fs, coloring.ConstructColorer(), labeler);
}

//fixme code duplication
template<class Graph>
void PrintColoredGraphWithColorFilter(const Graph &g, const ColorHandler<Graph> &coloring,
    const EdgesPositionHandler<Graph> &pos, const string &output_folder) {

  size_t edge_length_bound = 1000000;
  size_t colors_number = coloring.max_colors();
  TColorSet restricted_color = TColorSet::AllColorsSet(colors_number);

    shared_ptr<omnigraph::GraphSplitter<Graph>> splitter = ReliableSplitter<Graph>(g, edge_length_bound, 30);
    shared_ptr<omnigraph::GraphComponentFilter<Graph>> filter = make_shared<ComponentSingleColorFilter<Graph>>(g, coloring, restricted_color, edge_length_bound, 2);
    shared_ptr<omnigraph::GraphSplitter<Graph>> fs = make_shared<omnigraph::FilteringSplitterWrapper<Graph>>(splitter, filter);
    LengthIdGraphLabeler<Graph> basic_labeler(g);
    EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

    visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
    visualization::visualization_utils::WriteComponents(g, output_folder, fs, coloring.ConstructColorer(), labeler);
}

//todo alert!!! magic constants!!!
//todo refactoring of params needed
template<class gp_t>
void WriteComponentsAlongSequence(
        const gp_t& gp,
        const AbstractFilter<vector<typename gp_t::graph_t::VertexId>>& /*filter*/,
        const string& /*file_name*/,
        size_t /*split_edge_length*/, size_t /*component_vertex_number*/,
        const Sequence& /*s*/, const ColorHandler<typename gp_t::graph_t>& /*coloring*/) {
    typedef typename gp_t::graph_t Graph;
    LengthIdGraphLabeler < Graph > basic_labeler(gp.g);
    EdgePosGraphLabeler < Graph > pos_labeler(gp.g, gp.edge_pos);
    visualization::graph_labeler::CompositeLabeler < Graph > labeler(basic_labeler, pos_labeler);
}

template<class gp_t>
void PrintColoredGraphAlongRef(const gp_t& gp,
        const ColorHandler<Graph>& coloring,
        const string& output_filename) {
    LengthIdGraphLabeler < Graph > basic_labeler(gp.g);
    EdgePosGraphLabeler < Graph > pos_labeler(gp.g, gp.edge_pos);

    visualization::graph_labeler::CompositeLabeler < Graph > labeler(basic_labeler, pos_labeler);

//      only breakpoints
    TrivialBreakpointFinder<Graph> bp_f(gp.g, coloring, gp.edge_pos);

    WriteComponentsAlongSequence(gp, bp_f, labeler, output_filename, 1000000,
            30, MapperInstance(gp)->MapSequence(gp.genome.GetSequence()),
            *ConstructBorderColorer(gp.g, coloring)
//              *ConstructColorer(coloring)
                    );
}

}

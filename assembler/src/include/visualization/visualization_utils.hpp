#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "graph_printer.hpp"
#include "omni/omni_utils.hpp"
#include "dijkstra/dijkstra_helper.hpp"
#include "omni/splitters.hpp"
#include "omni/graph_component.hpp"
#include "visualizers.hpp"
#include "vertex_linker.hpp"

namespace omnigraph {
namespace visualization {


template<class Graph>
void WriteComponents(const Graph& g,
        const string& folder_name,
        shared_ptr<GraphSplitter<Graph>> inner_splitter,
        shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
//  shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
    auto filter = make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CollectingSplitterWrapper<Graph>>(inner_splitter, filter);
    omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(*splitter, folder_name);
}

template<class Graph>
void DrawComponentsOfShortEdges(const Graph& g, size_t min_length, size_t sinks, size_t sources)
{
    vector<typename Graph::EdgeId> short_edges;
    std::string pics_folder_ = cfg::get().output_dir + ToString(min_length) + "_" + ToString(sinks) + "_" + ToString(sources) + "_"+ "pics_polymorphic/";
    make_dir(pics_folder_);
    INFO("Writing pics with components consisting of short edges to " + pics_folder_);
    shared_ptr<GraphSplitter<Graph>> splitter = LongEdgesExclusiveSplitter<Graph>(g, min_length);
    while (splitter->HasNext()) {
        GraphComponent<Graph> component = splitter->Next();
        if(component.v_size() > 3 && component.sinks().size() == sinks && component.sources().size() == sources)
        {
        	bool fail = false;
        	for(auto v : component.sources()) {
        		if(component.g().IncomingEdgeCount(v) != 1) {
        			fail = true;
        		}
        	}
        	for(auto v : component.sinks()) {
        		if(component.g().OutgoingEdgeCount(v) != 1) {
        			fail = true;
        		}
        	}

        	if(fail)
        	{
        		continue;
        	}

            StrGraphLabeler<Graph> labeler(component.g());
            CoverageGraphLabeler<Graph> labeler2(component.g());
            CompositeLabeler<Graph> compositeLabeler(labeler, labeler2);
            WriteComponentSinksSources(component, pics_folder_ + ToString(g.int_id(*component.vertices().begin()))
                                                                                   + ".dot", visualization::DefaultColorer(g),
                                                                                   compositeLabeler);
            INFO("Component is written to " + ToString(g.int_id(*component.vertices().begin())) + ".dot");

            //            PrintComponent(component,
//                                pics_folder_ + "ShortComponents/"
//                                        + ToString(gp.g.int_id(component.vertices_[0]))
//                                         + ".dot");
        }
    }
}


template<class Graph>
void WriteSizeLimitedComponents(const Graph& g,
        const string& folder_name,
        shared_ptr<GraphSplitter<Graph>> inner_splitter,
        shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler, int min_component_size, int max_component_size, size_t max_components) {
    EmptyGraphLinker<Graph> linker;

    auto filter = make_shared<omnigraph::ComponentSizeFilter<Graph>>(g, 1000000000, (size_t) min_component_size, (size_t) max_component_size);
    shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CollectingSplitterWrapper<Graph>>(inner_splitter, filter);
    omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker, false, max_components).SplitAndVisualize(*splitter, folder_name);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph>& gc,
        const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    BorderDecorator<Graph> component_colorer(gc, *colorer, "yellow");
    ofstream os;
    os.open(file_name);
    omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), true).Visualize(gc, os, labeler, component_colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentSinksSources(const GraphComponent<Graph>& gc,
        const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    SinkSourceDecorator<Graph> component_colorer(gc, *colorer);
    ofstream os;
    os.open(file_name);
    omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), true).Visualize(gc, os, labeler, component_colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentSinksSources(const GraphComponent<Graph>& gc,
        const string& file_name) {

    StrGraphLabeler<Graph> labeler(gc.g());
    CoverageGraphLabeler<Graph> labeler2(gc.g());
    CompositeLabeler<Graph> compositeLabeler(labeler, labeler2);
    EmptyGraphLinker<Graph> linker;
    WriteComponentSinksSources(gc, file_name, DefaultColorer(gc.g()),
            compositeLabeler);
}

template<class Graph>
void WriteSimpleComponent(const GraphComponent<Graph>& gc,
        const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    ofstream os;
    os.open(file_name);
    omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), false).Visualize(gc, os, labeler, *colorer, linker);
    os.close();
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g, vector<typename Graph::EdgeId> path,
        const string& prefix_path, shared_ptr<GraphColorer<Graph>> colorer,
        const GraphLabeler<Graph> &labeler, bool color_path = true) {
    auto edge_colorer = make_shared<CompositeEdgeColorer<Graph>>("black");
    edge_colorer->AddColorer(colorer);
    if (color_path) {
        edge_colorer->AddColorer(make_shared<SetColorer<Graph>>(g, path, "green"));
    }
    shared_ptr<GraphColorer<Graph>> resulting_colorer =  make_shared<CompositeGraphColorer<Graph>>(colorer, edge_colorer);
    shared_ptr<GraphSplitter<Graph>> rs = ReliableSplitterAlongPath<Graph>(g, path);
    auto filter = make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CondensingSplitterWrapper<Graph>>(rs, filter);
    WriteComponents<Graph>(g, prefix_path, splitter, resulting_colorer, labeler);
}

template<class Graph>
class LocalityPrintingRH {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& g_;
    const GraphLabeler<Graph>& labeler_;
    std::shared_ptr<visualization::GraphColorer<Graph>> colorer_;
    const string output_folder_;
public:
    LocalityPrintingRH(const Graph& g
            , const GraphLabeler<Graph>& labeler
            , std::shared_ptr<visualization::GraphColorer<Graph>> colorer
            , const string& output_folder) :
            g_(g),
            labeler_(labeler),
            colorer_(colorer),
            output_folder_(output_folder) {
        path::make_dirs(output_folder_);
    }

    void HandleDelete(EdgeId e, const string& add_label = "") {
        //todo magic constant
//          map<EdgeId, string> empty_coloring;
        auto edge_colorer = make_shared<visualization::CompositeEdgeColorer<Graph>>("black");
        edge_colorer->AddColorer(colorer_);
        edge_colorer->AddColorer(make_shared<visualization::SetColorer<Graph>>(this->g(), vector<EdgeId>(1, e), "green"));
        shared_ptr<visualization::GraphColorer<Graph>> resulting_colorer = make_shared<visualization::CompositeGraphColorer<Graph>>(colorer_, edge_colorer);

        string fn = output_folder_ + "edge_" + ToString(this->g().int_id(e)) + add_label + ".dot";
        omnigraph::visualization::WriteComponent(omnigraph::EdgeNeighborhood<Graph>(this->g(), e, 50, 250)
                , fn
                , resulting_colorer, labeler_);
    }

private:
    DECL_LOGGER("LocalityPrintingRH")
    ;
};

//static void WriteFilteredComponents(const Graph& g,
//      const string& folder_name,
//      shared_ptr<GraphComponentFilter<Graph>> filter,
//      shared_ptr<GraphSplitter<Graph>> splitter,
//      shared_ptr<GraphColorer<Graph>> colorer,
//      const GraphLabeler<Graph> &labeler) {
//  EmptyGraphLinker<Graph> linker;
////    shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
//  omnigraph::FilteringSplitterWrapper<Graph> filtered_splitter(splitter, filter);
//  omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(filtered_splitter, folder_name);
//}

}
}

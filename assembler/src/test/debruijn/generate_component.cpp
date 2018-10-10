//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "io/binary/graph_pack.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

using namespace std;

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

template<class Graph>
class BlockedComponentFinder : public AbstractNeighbourhoodFinder<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    set<VertexId> blocking_vertices_;
    const size_t edge_length_bound_;

    VertexId OtherEnd(EdgeId e, VertexId v) const {
        if (this->graph().EdgeStart(e) == v)
            return this->graph().EdgeEnd(e);
        else
            return this->graph().EdgeStart(e);
    }

    void Go(VertexId v, set<VertexId>& entered) const {
        if (blocking_vertices_.count(v) || entered.count(v))
            return;

        entered.insert(v);

        for (EdgeId e : this->graph().IncidentEdges(v)) {
            VertexId adjacent_v = OtherEnd(e, v);
            if (this->graph().length(e) <= edge_length_bound_) {
                Go(adjacent_v, entered);
            }
        }
    }

public:

    BlockedComponentFinder(const Graph &graph, const vector<VertexId>& blocking_vertices, size_t edge_length_bound)
            : AbstractNeighbourhoodFinder<Graph>(graph),
              blocking_vertices_(blocking_vertices.begin(), blocking_vertices.end()),
              edge_length_bound_(edge_length_bound) {
    }

    GraphComponent<Graph> Find(typename Graph::VertexId v) const {
        set<VertexId> result;
        Go(v, result);
        return GraphComponent<Graph>::FromVertices(this->graph(), result);
    }

    vector<VertexId> InnerVertices(const GraphComponent<Graph>& /*component*/) const {
        VERIFY(false);
        return vector<VertexId>();
    }
};

void Launch(size_t K, string saves_path, size_t start_vertex_int_id,
            vector<size_t> blocking_int_ids, size_t edge_length_bound,
            string component_out_path) {
    conj_graph_pack gp(K, "tmp", 0);
    omnigraph::GraphElementFinder<Graph> element_finder(gp.g);
    io::binary::BasePackIO<Graph> io;
    io.Load(saves_path, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");
    VertexId starting_vertex = element_finder.ReturnVertexId(start_vertex_int_id);
    vector<VertexId> blocking_vertices;
    for (size_t int_id : blocking_int_ids) {
        blocking_vertices.push_back(element_finder.ReturnVertexId(int_id));
    }
    BlockedComponentFinder<Graph> component_finder(gp.g, blocking_vertices, edge_length_bound);
    GraphComponent<Graph> component_to_save = ComponentCloser<Graph>(gp.g, 0).CloseComponent(component_finder.Find(starting_vertex));
    INFO("Blocked component has " << component_to_save.v_size() << " vertices");
    io.Save(component_out_path, gp);
    gp.edge_pos.Attach();
    visualization::visualization_utils::WriteComponent<Graph>(component_to_save, component_out_path + ".dot", debruijn_graph::stats::DefaultColorer(gp),
                                                    visualization::graph_labeler::DefaultLabeler<Graph>(gp.g, gp.edge_pos));
}

}

int main(int argc, char** argv) {
    if (argc < 6) {
        cout << "Usage: component_generator <K> "
             << "<saves path> <start vertex id> <edge length bound> "
             << "<component output> [<blocking vertex ids>]" << endl;
        exit(1);
    }
//    TmpFolderFixture("tmp");
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    size_t start_vertex_int_id = std::stoll(argv[3]);
    INFO("Start vertex " << start_vertex_int_id);
    size_t edge_length_bound = std::stoll(argv[4]);
    INFO("Edge length bound " << edge_length_bound);
    string component_out_path = argv[5];
    INFO("Save component to " << component_out_path);
    vector<size_t> blocking_int_ids;
    for (int i = 6; i < argc; ++i) {
        blocking_int_ids.push_back(std::stoll(argv[i]));
    }
    INFO("Blocking ids " << blocking_int_ids);
    debruijn_graph::Launch(K, saves_path, start_vertex_int_id, blocking_int_ids, edge_length_bound, component_out_path);
    return 0;
}

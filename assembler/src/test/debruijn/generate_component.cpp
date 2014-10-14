#include "standard_base.hpp"
#include "logger/log_writers.hpp"

#include "graphio.hpp"
#include "graph_pack.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {
    const string tmp_folder = "tmp/";

    template<class Graph>
    class BlockedComponentFinder : public AbstractNeighbourhoodFinder<Graph> {
    private:
        typedef typename Graph::VertexId VertexId;
        typedef typename Graph::EdgeId EdgeId;

        set<VertexId> blocking_vertices_;
        const size_t edge_length_bound_;

        VertexId OtherEnd(EdgeId e, VertexId v) {
            if (this->graph().EdgeStart(e) == v)
                return this->graph().EdgeEnd(e);
            else
                return this->graph().EdgeStart(e);
        }

        void Go(VertexId v, set<VertexId>& entered) {
            if (blocking_vertices_.count(v) || entered.count(v))
                return;

            entered.insert(v);

            for (EdgeId e : this->graph().AdjacentEdges(v)) {
                VertexId adjacent_v = OtherEnd(e, v);
                if (this->graph().length(e) <= edge_length_bound_) {
                    Go(adjacent_v, entered);
                }
            }
        }

    public:

        BlockedComponentFinder(const Graph &graph, const vector<VertexId>& blocking_vertices,
                               size_t edge_length_bound)
                : AbstractNeighbourhoodFinder<Graph>(graph),
                  blocking_vertices_(blocking_vertices.begin(), blocking_vertices.end()),
                  edge_length_bound_(edge_length_bound) {
        }

        GraphComponent<Graph> Find(typename Graph::VertexId v) {
            set<VertexId> result;
            Go(v, result);
            return GraphComponent<Graph>(this->graph(), result.begin(), result.end());
        }

        vector<VertexId> InnerVertices(const GraphComponent<Graph>& /*component*/) {
            VERIFY(false);
            return vector<VertexId>();
        }
    };

    void Launch(size_t K, string saves_path,
                size_t start_vertex_int_id, vector<size_t> blocking_int_ids,
                size_t edge_length_bound, string component_out_path) {
        conj_graph_pack gp(55, tmp_folder, 0);
        graphio::ScanGraphPack(saves_path, gp);
        omnigraph::GraphElementFinder<Graph> element_finder(gp.g);
        VertexId starting_vertex = element_finder.ReturnVertexId(start_vertex_int_id);
        vector<VertexId> blocking_vertices;
        for (size_t int_id : blocking_int_ids) {
            blocking_vertices.push_back(element_finder.ReturnVertexId(int_id));
        }
        BlockedComponentFinder<Graph> component_finder(gp.g, blocking_vertices, edge_length_bound);
        GraphComponent<Graph> component_to_save = component_finder.Find(starting_vertex);
        graphio::ConjugateDataPrinter<Graph> printer(component_to_save);
        graphio::PrintGraphPack(component_out_path, printer, gp);
    }
}

int main(int /*argc*/, char** argv) {
    create_console_logger();
    size_t K = 55;
    string saves_path = "";
    size_t start_vertex_int_id = 0;
    vector<size_t> blocking_int_ids = {};
    size_t edge_length_bound = 1000;
    string component_out_path;

    debruijn_graph::Launch(K, saves_path, start_vertex_int_id, blocking_int_ids, edge_length_bound, component_out_path);
    return 0;
}

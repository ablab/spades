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
#include "assembly_graph/components/splitters.hpp"

using namespace std;

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

std::string EdgeName(const Graph& g, EdgeId e) {
    return io::MakeContigId(g.int_id(e), g.length(e) + g.k(), g.coverage(e), "EDGE");
}

std::string CanonicalName(const Graph& g, EdgeId e) {
    if (e <= g.conjugate(e)) {
        return EdgeName(g, e);
    } else {
        return EdgeName(g, g.conjugate(e)) + "'";
    }
}

template<class EdgeContainer>
std::string NeighboursString(const Graph& g, const EdgeContainer& edges) {
    std::stringstream ss;
    string delim = "";
    for (EdgeId e: edges) {
        ss << delim << CanonicalName(g, e);
        delim = ",";
    }
    return ss.str();
}

void Launch(size_t K, string saves_path, size_t edge_length_bound,
            const string& fastg_output) {
    fs::TmpFolderFixture tmp_dir("tmp");
    //TODO no need for whole graph pack; change to Graph
    conj_graph_pack gp(K, "tmp", 0);
    io::binary::BasePackIO<Graph>().Load(saves_path, gp);

    io::OFastaReadStream oss(fastg_output);
    for (auto it = gp.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        if (gp.g.length(e) > edge_length_bound) {
            DEBUG("Processing edge " << gp.g.str(e));

            typedef ComposedDijkstraSettings<Graph,
                    BoundedEdgeLenCalculator<Graph>,
                    ZeroLengthProcessChecker<Graph>,
                    VertexPutChecker<Graph>,
                    ForwardNeighbourIteratorFactory<Graph>> ForwardShortEdgeDijkstraSettings;

            typedef Dijkstra<Graph, ForwardShortEdgeDijkstraSettings> ForwardShortEdgeDijkstra;

            ForwardShortEdgeDijkstraSettings settings(
                    BoundedEdgeLenCalculator<Graph>(gp.g, edge_length_bound),
                    ZeroLengthProcessChecker<Graph>(),
                    VertexPutChecker<Graph>(),
                    ForwardNeighbourIteratorFactory<Graph>(gp.g));
            ForwardShortEdgeDijkstra dijkstra(gp.g, settings,
                                              std::numeric_limits<size_t>::max());
            dijkstra.Run(gp.g.EdgeEnd(e));

            set<EdgeId> long_reachable;
            for (VertexId v : dijkstra.ProcessedVertices()) {
                for (EdgeId out_e : gp.g.OutgoingEdges(v)) {
                    if (gp.g.length(out_e) > edge_length_bound) {
                        long_reachable.insert(out_e);
                    }
                }
            }

            for (EdgeId n : long_reachable) {
                DEBUG("Neighbour " << gp.g.str(n));
                auto path = dijkstra.GetShortestPathTo(gp.g.EdgeStart(n));
                DEBUG("Connecting path " << PrintPath(gp.g, path));
            }

            std::string neighbours_string = NeighboursString(gp.g, long_reachable);

            string name = CanonicalName(gp.g, e);
            if (!neighbours_string.empty()) {
                name += ":" + neighbours_string;
            }
            name += ";";
            oss << io::SingleRead(name, gp.g.EdgeNucls(e).str());
        }
    }
}

//void Launch(size_t K, string saves_path, size_t edge_length_bound,
//            const string& fastg_output) {
//    TmpFolderFixture tmp_dir("tmp");
//    //TODO no need for whole graph pack; change to Graph
//    conj_graph_pack gp(K, "tmp", 0);
//    graphio::ScanGraphPack(saves_path, gp);
//    auto splitter = omnigraph::LongEdgesExclusiveSplitter(gp.g, edge_length_bound);
//    io::osequencestream oss(fastg_output);
//    while (splitter->HasNext()) {
//        auto gc = splitter->Next();
//        set<EdgeId> long_outgoing;
//        set<EdgeId> long_incoming;
//        for (auto v : gc.vertices()) {
//            for (auto e : gp.g.OutgoingEdges(v))
//                if (gp.g.length(e) > edge_length_bound)
//                    long_outgoing.insert(e);
//            for (auto e : gp.g.IncomingEdges(v))
//                if (gp.g.length(e) > edge_length_bound)
//                    long_incoming.insert(e);
//        }
//        std::stringstream ss;
//        string delim = "";
//        for (EdgeId e: long_outgoing) {
//            ss << delim << CanonicalName(gp.g, e);
//            delim = ",";
//        }
//        string outgoing_names = ss.str();
//        for (EdgeId e : long_incoming) {
//            string name = CanonicalName(gp.g, e);
//            if (!outgoing_names.empty()) {
//                name += ":" + outgoing_names;
//            }
//            name += ";";
//            oss << io::SingleRead(name, gp.g.EdgeNucls(e).str());
//        }
//    }
//}

}

int main(int argc, char** argv) {
    if (argc < 5) {
        cout << "Usage: short_edge_contractor <K> "
             << "<saves path> <edge length bound> "
             << "<fastg output>" << endl;
        exit(1);
    }
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    size_t edge_length_bound = std::stoll(argv[3]);
    INFO("Edge length bound " << edge_length_bound);
    string fastg_output = argv[4];
    debruijn_graph::Launch(K, saves_path, edge_length_bound, fastg_output);
    return 0;
}

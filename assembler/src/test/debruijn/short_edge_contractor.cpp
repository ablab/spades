//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "utils/simple_tools.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "assembly_graph/components/splitters.hpp"

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

void Launch(size_t K, string saves_path, size_t edge_length_bound,
            const string& fastg_output) {
    TmpFolderFixture tmp_dir("tmp");
    //TODO no need for whole graph pack; change to Graph
    conj_graph_pack gp(K, "tmp", 0);
    graphio::ScanGraphPack(saves_path, gp);
    auto splitter = omnigraph::LongEdgesExclusiveSplitter(gp.g, edge_length_bound);
    io::osequencestream oss(fastg_output);
    while (splitter->HasNext()) {
        auto gc = splitter->Next();
        set<EdgeId> long_outgoing;
        set<EdgeId> long_incoming;
        for (auto v : gc.vertices()) {
            for (auto e : gp.g.OutgoingEdges(v))
                if (gp.g.length(e) > edge_length_bound)
                    long_outgoing.insert(e);
            for (auto e : gp.g.IncomingEdges(v))
                if (gp.g.length(e) > edge_length_bound)
                    long_incoming.insert(e);
        }
        std::stringstream ss;
        string delim = "";
        for (EdgeId e: long_outgoing) {
            ss << delim << CanonicalName(gp.g, e);
            delim = ",";
        }
        string outgoing_names = ss.str();
        for (EdgeId e : long_incoming) {
            string name = CanonicalName(gp.g, e);
            if (!outgoing_names.empty()) {
                name += ":" + outgoing_names;
            }
            name += ";";
            oss << io::SingleRead(name, gp.g.EdgeNucls(e).str());
        }
    }
}

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

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "standard_base.hpp"
#include "simple_tools.hpp"
#include "logger/log_writers.hpp"

#include "graphio.hpp"
#include "graph_pack.hpp"
#include "io/io_helper.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

void Launch(size_t k, const string& saves_path, const string& contigs_path) {
    conj_graph_pack gp(k, "tmp", 0);
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);

    auto contigs_stream_ptr = io::EasyStream(contigs_path, false);
    auto mapper_ptr = MapperInstance(gp);

    io::SingleRead contig;
    while (!contigs_stream_ptr->eof()) {
        (*contigs_stream_ptr) >> contig;

        string name = contig.name();
        vector<EdgeId> path = mapper_ptr->MapRead(contig).simple_path();

        /*Your code here*/
    }
}

}

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "Usage: annotation_propagator <K> <saves path> <contigs_path>"  << endl;
        exit(1);
    }
    TmpFolderFixture("tmp");
    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string saves_path = argv[2];
    string contigs_path = argv[3];
    debruijn_graph::Launch(k, saves_path, contigs_path);
    return 0;
}

//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/simple_tools.hpp"
#include "dev_support/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "io/reads_io/file_reader.hpp"
#include "propagate.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

int main(int argc, char** argv) {
    using namespace debruijn_graph;

    if (argc < 6) {
        cout << "Usage: annotation_propagator <K> <config path> <saves path> <contigs_path> "
                "<init binning info> <final binning info> (<bins of interest>)*"  << endl;
        exit(1);
    }
    //TmpFolderFixture fixture("tmp");
    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string config_path = argv[2];
    string saves_path = argv[3];
    string contigs_path = argv[4];
    string annotation_in_fn = argv[5];
    string annotation_out_fn = argv[6];
//    debruijn_graph::Launch(k, saves_path, contigs_path);

    std::vector<bin_id> bins_of_interest;
    for (int i = 7; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    cfg::create_instance(config_path);
    conj_graph_pack gp(k, "tmp", cfg::get().ds.reads.lib_count());
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanWithClusteredIndices(saves_path, gp, gp.clustered_indices);
    auto contigs_stream_ptr = make_shared<io::FileReadStream>(contigs_path);

    AnnotationPropagator propagator(gp);
    propagator.Run(*contigs_stream_ptr, annotation_in_fn, bins_of_interest, annotation_out_fn);

    return 0;
}

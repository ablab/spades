//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/logger/log_writers.hpp"
#include "pipeline/graphio.hpp"
#include "read_binning.hpp"
#include "propagate.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

int main(int argc, char** argv) {
    using namespace debruijn_graph;

    if (argc < 9) {
        cout << "Usage: prop_binning <K> <config path> <saves path> <contigs path> <binning info> "
                "<left reads> <right reads> <output root> <sample name> (<bins of interest>)*"  << endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string config_path = argv[2];
    string saves_path = argv[3];
    string contigs_path = argv[4];
    string annotation_path = argv[5];
    //TODO: don't save the propagated info
    string contigs_binning_path = annotation_path + ".prop";
    string left_reads = argv[6];
    string right_reads = argv[7];
    string out_root = argv[8];
    string sample_name = argv[9];

    std::vector<bin_id> bins_of_interest;
    for (int i = 10; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    cfg::create_instance(config_path);
    conj_graph_pack gp(k, "tmp", cfg::get().ds.reads.lib_count());
    gp.kmer_mapper.Attach();

    INFO("Load graph from " << saves_path);
    graphio::ScanWithClusteredIndices(saves_path, gp, gp.clustered_indices);

    //Propagation stage
    io::FileReadStream contigs_stream(contigs_path);
    AnnotationPropagator propagator(gp);
    propagator.Run(contigs_stream, annotation_path, bins_of_interest, contigs_binning_path);

    //Binning stage
    ContigBinner binner(gp, bins_of_interest);
    AnnotationStream binning_stream(contigs_binning_path);
    contigs_stream.reset();
    binner.Init(out_root, sample_name, contigs_stream, binning_stream);

    auto paired_stream = io::PairedEasyStream(left_reads, right_reads, false, 0);
    binner.Run(*paired_stream);
    binner.close();

    return 0;
}

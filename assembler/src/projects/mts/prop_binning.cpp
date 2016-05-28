//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "getopt_pp/getopt_pp.h"
#include "dev_support/logger/log_writers.hpp"
#include "pipeline/graphio.hpp"
#include "read_binning.hpp"
#include "propagate.hpp"
#include <modules/visualization/position_filler.hpp>

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

int main(int argc, char** argv) {
    using namespace debruijn_graph;
    using namespace GetOpt;

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    size_t k;
    string saves_path, contigs_path, annotation_path;
    string left_reads, right_reads;
    string out_root, sample_name;
    std::vector<bin_id> bins_of_interest;
    bool no_binning;
    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('c', contigs_path)
            >> Option('a', annotation_path)
            >> Option('l', left_reads)
            >> Option('r', right_reads)
            >> Option('o', out_root)
            >> Option('n', sample_name)
            >> Option('b', bins_of_interest, {})
            >> OptionPresent('p', no_binning);
    } catch(GetOptEx &ex) {
        cout << "Usage: prop_binning -k <K> -s <saves path> -c <contigs path> -a <binning annotation> "
                "-l <left reads> -r <right reads> -o <output root> -n <sample name> "
                "[-p to disable binning] [-b <bins of interest>*]"  << endl;
        exit(1);
    }
    //TODO: don't save the propagated info
    string propagated_path = annotation_path + ".prop";

    conj_graph_pack gp(k, "tmp", 1);
    gp.kmer_mapper.Attach();

    INFO("Load graph and clustered paired info from " << saves_path);
    graphio::ScanWithClusteredIndices(saves_path, gp, gp.clustered_indices);

//    gp.edge_pos.Attach();
//    FillPos(gp, "/Sid/snurk/mts/data/infant_gut/refs/Enterococcus_faecalis_37_1.fasta",
//            "CAG1_ref_", true);

    //Propagation stage
    //TODO: make this optional
    INFO("Using contigs from " << contigs_path);
    io::FileReadStream contigs_stream(contigs_path);
    INFO("Propagation launched");
    AnnotationPropagator propagator(gp);
    propagator.Run(contigs_stream, annotation_path, bins_of_interest, propagated_path);
    INFO("Propagation finished");

    if (no_binning) {
        INFO("Binning was disabled with -p flag");
        return 0;
    }
    //Binning stage
    contigs_stream.reset();
    ContigBinner binner(gp, bins_of_interest);
    INFO("Initializing binner");
    INFO("Using propagated annotation from " << propagated_path);
    AnnotationStream binning_stream(propagated_path);
    binner.Init(out_root, sample_name, contigs_stream, binning_stream);

    auto paired_stream = io::PairedEasyStream(left_reads, right_reads, false, 0);
    INFO("Running binner on " << left_reads << " and " << right_reads);
    binner.Run(*paired_stream);
    binner.close();

    return 0;
}

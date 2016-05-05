//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/simple_tools.hpp"
#include "dev_support/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "io/reads_io/file_reader.hpp"
#include "read_binning.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

void ContigBinner::Init(const string& output_root, const string& sample_name, io::SingleStream& contigs, AnnotationStream& annotation_stream) {
    edge_annotation_.Fill(contigs, annotation_stream);
    for (bin_id bin : edge_annotation_.interesting_bins()) {
        string out_dir = output_root + "/" + ToString(bin) + "/";
        path::make_dirs(out_dir);
        out_streams_.insert(make_pair(bin,
                                      make_shared<io::OPairedReadStream>(out_dir + sample_name + "_1.fastq",
                                                                                  out_dir + sample_name + "_2.fastq")));
    }
}

void ContigBinner::Run(io::PairedStream& paired_reads) {
    io::PairedRead paired_read;
    while (!paired_reads.eof()) {
        paired_reads >> paired_read;
        set<bin_id> bins;
        insert_all(bins, edge_annotation_.RelevantBins(paired_read.first()));
        insert_all(bins, edge_annotation_.RelevantBins(paired_read.second()));
        for (auto bin : bins) {
            (*(out_streams_[bin])) << paired_read;
        }
    }
}

};

//todo make it take dataset info
int main(int argc, char** argv) {
    using namespace debruijn_graph;

    if (argc < 9) {
        cout << "Usage: read_binning <K> <saves path> <contigs path> <contigs binning info> "
                "<left reads> <right reads> <output root> <sample name> (<bins of interest>)*"  << endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string saves_path = argv[2];
    string contigs_path = argv[3];
    string contigs_binning_path = argv[4];
    string left_reads = argv[5];
    string right_reads = argv[6];
    string out_root = argv[7];
    string sample_name = argv[8];

    std::vector<bin_id> bins_of_interest;
    for (int i = 9; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    conj_graph_pack gp(k, "tmp", 0);
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);

    ContigBinner binner(gp, bins_of_interest);

    auto contigs_stream_ptr = make_shared<io::FileReadStream>(contigs_path);
    AnnotationStream binning_stream(contigs_binning_path);

    binner.Init(out_root, sample_name, *contigs_stream_ptr, binning_stream);

    auto paired_stream = io::PairedEasyStream(left_reads, right_reads, false, 0);
    binner.Run(*paired_stream);
    binner.close();
    return 0;
}

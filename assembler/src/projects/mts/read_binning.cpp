//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/simple_tools.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "io/reads/file_reader.hpp"
#include "read_binning.hpp"

namespace debruijn_graph {

set<bin_id> ContigBinner::RelevantBins(const io::SingleRead& r) const {
    return edge_annotation_.RelevantBins(mapper_->MapRead(r).simple_path());
}

void ContigBinner::Init(bin_id bin) {
    string out_dir = out_root_ + "/" + ToString(bin) + "/";
    path::make_dirs(out_dir);
    out_streams_.insert(make_pair(bin, make_shared<io::OPairedReadStream>(out_dir + sample_name_ + "_1.fastq",
                                                                          out_dir + sample_name_ + "_2.fastq")));
}

void ContigBinner::Run(io::PairedStream& paired_reads) {
    io::PairedRead paired_read;
    while (!paired_reads.eof()) {
        paired_reads >> paired_read;
        set<bin_id> bins;
        insert_all(bins, RelevantBins(paired_read.first()));
        insert_all(bins, RelevantBins(paired_read.second()));
        for (auto bin : bins) {
            if (out_streams_.find(bin) == out_streams_.end()) {
                Init(bin);
            }
            (*(out_streams_[bin])) << paired_read;
        }
    }
}

};

//todo make it take dataset info
/*
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
*/

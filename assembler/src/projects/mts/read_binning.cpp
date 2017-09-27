//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "io/reads/file_reader.hpp"
#include "read_binning.hpp"

namespace debruijn_graph {

set<bin_id> ContigBinner::RelevantBins(const io::SingleRead& r) const {
    return edge_annotation_.RelevantBins(mapper_->MapRead(r).simple_path());
}

void ContigBinner::Init(bin_id bin) {
    string out_dir = out_root_ + "/" + bin + "/";
    fs::make_dirs(out_dir);
    out_streams_.insert(make_pair(bin, make_shared<ContigBinner::Stream>(
        out_dir + sample_name_ + "_1.fastq.gz",
        out_dir + sample_name_ + "_2.fastq.gz")
    ));
}

void ContigBinner::Run(io::PairedStream& paired_reads) {
    io::PairedRead paired_read;
    while (!paired_reads.eof()) {
        paired_reads >> paired_read;
        set<bin_id> bins;
        utils::insert_all(bins, RelevantBins(paired_read.first()));
        utils::insert_all(bins, RelevantBins(paired_read.second()));
        for (const auto& bin : bins) {
            if (bins_of_interest_.size() && !bins_of_interest_.count(bin)) {
                INFO(bin << " was excluded from read binning");
                continue;
            }
            if (out_streams_.find(bin) == out_streams_.end()) {
                Init(bin);
            }
            (*(out_streams_[bin])) << paired_read;
        }
    }
}

void BinReads(const conj_graph_pack& gp, const std::string& out_root,
             const std::string& sample,
             const std::string& left_reads, const std::string& right_reads,
             const EdgeAnnotation& edge_annotation,
             const vector<string>& bins_of_interest) {
    ContigBinner binner(gp, edge_annotation, out_root, sample, bins_of_interest);
    INFO("Initializing binner for " << sample);
    auto paired_stream = io::PairedEasyStream(left_reads, right_reads, false, 0);
    INFO("Running binner on " << left_reads << " and " << right_reads);
    binner.Run(*paired_stream);
}

};

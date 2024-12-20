//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "io/reads/file_reader.hpp"
#include "read_binning.hpp"

namespace debruijn_graph {

std::set<bin_id> ContigBinner::RelevantBins(const io::SingleRead& r) const {
    return edge_annotation_.RelevantBins(mapper_->MapRead(r).simple_path());
}

void ContigBinner::Init(bin_id bin) {
    std::filesystem::path out_dir = out_root_ / bin;
    create_directories(out_dir);
    out_streams_.emplace(bin, std::make_unique<ContigBinner::Stream>(
        out_dir / sample_name_.concat("_1.fastq.gz"),
        out_dir / sample_name_.concat("_2.fastq.gz"))
    );
}

void ContigBinner::Run(io::PairedStream& paired_reads) {
    io::PairedRead paired_read;
    while (!paired_reads.eof()) {
        paired_reads >> paired_read;
        std::set<bin_id> bins;
        utils::insert_all(bins, RelevantBins(paired_read.first()));
        utils::insert_all(bins, RelevantBins(paired_read.second()));
        for (const auto& bin : bins) {
            if (bins_of_interest_.size() && !bins_of_interest_.count(bin)) {
                INFO(bin << " was excluded from read binning");
                continue;
            }
            if (!out_streams_.count(bin)) {
                Init(bin);
            }
            *out_streams_[bin] << paired_read;
        }
    }
}

void BinReads(const graph_pack::GraphPack& gp, const std::string& out_root,
              const std::string& sample,
              const std::string& left_reads, const std::string& right_reads,
              const EdgeAnnotation& edge_annotation,
              const std::vector<std::string>& bins_of_interest) {
    ContigBinner binner(gp, edge_annotation, out_root, sample, bins_of_interest);
    INFO("Initializing binner for " << sample);
    auto paired_stream = io::PairedEasyStream(left_reads, right_reads, false, 0);
    INFO("Running binner on " << left_reads << " and " << right_reads);
    binner.Run(paired_stream);
}

};

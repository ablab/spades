//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_splitting.hpp"
#include "binning.hpp"

#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "modules/alignment/kmer_sequence_mapper.hpp"

#include "io/reads/osequencestream.hpp"
#include "io/reads/io_helper.hpp"
#include "io/dataset_support/read_converter.hpp"

#include <threadpool/threadpool.hpp>

#include <vector>

using namespace debruijn_graph;
using namespace bin_stats;

namespace binning {

class ReadsPathsToComponentsListener : public debruijn_graph::SequenceMapperListener {
 public:
    ReadsPathsToComponentsListener(const Binning &binning,
                                   std::vector<io::OFastqPairedStream> &binned_reads_ostreams,
                                   io::OFastqPairedStream &unbinned_reads_ostream,
                                   const SoftBinsAssignment &soft_edge_labels,
                                   const BinningAssignmentStrategy &assignment_strategy,
                                   const double bin_weight_threshold)
        : binned_reads_ostreams_(binned_reads_ostreams),
          unbinned_reads_ostream_(unbinned_reads_ostream),
          binning_(binning),
          soft_edge_labels_(soft_edge_labels),
          assignment_strategy_(assignment_strategy),
          bin_weight_threshold_(bin_weight_threshold) {}

    void StartProcessLibrary(size_t threads_count) override {
        binned_reads_.resize(threads_count, BinnedReads(binned_reads_ostreams_.size()));
        unbinned_reads_.resize(threads_count);
    }
    
    void MergeBuffer(size_t thread_index) override {
        BinnedReads &reads = binned_reads_[thread_index];
        ReadVector &unbinned =  unbinned_reads_[thread_index];

        for (size_t part_id = 0; part_id < binned_reads_ostreams_.size(); ++part_id) {
            io::OFastqPairedStream& paired_ostream = binned_reads_ostreams_[part_id];

            for (const auto& paired_read : reads[part_id])
                paired_ostream << paired_read;

            reads[part_id].clear();
        }

        for (const auto& paired_read : unbinned)
            unbinned_reads_ostream_ << paired_read;
        unbinned.clear();
    }

    void ProcessPairedRead(size_t thread_index,
                           const io::PairedReadSeq& paired_read,
                           const MappingPath<EdgeId>& first_path,
                           const MappingPath<EdgeId>& second_path) override {
        BinnedReads &reads = binned_reads_[thread_index];
        ReadVector &unbinned =  unbinned_reads_[thread_index];

        if (!first_path.empty() && !second_path.empty()) {
            auto first_read_bins =
                    assignment_strategy_.ChooseMajorBins(assignment_strategy_.AssignScaffoldBins(first_path.simple_path(),
                                                                                                 soft_edge_labels_, binning_),
                                                         soft_edge_labels_, binning_);
            auto second_read_bins =
                    assignment_strategy_.ChooseMajorBins(assignment_strategy_.AssignScaffoldBins(second_path.simple_path(),
                                                                                                 soft_edge_labels_, binning_),
                                                         soft_edge_labels_, binning_);

            std::unordered_set<bin_stats::Binning::BinId> read_bins;
            read_bins.insert(first_read_bins.begin(), first_read_bins.end());
            read_bins.insert(second_read_bins.begin(), second_read_bins.end());

            for (auto bin : read_bins)
                reads[bin].push_back(paired_read);
        } else 
            unbinned.push_back(paired_read);
    }

 private:
    using ReadVector = std::vector<io::PairedReadSeq>;
    using BinnedReads = std::vector<ReadVector>;
    // Per thread
    std::vector<BinnedReads> binned_reads_;
    std::vector<ReadVector> unbinned_reads_;

    std::vector<io::OFastqPairedStream> &binned_reads_ostreams_;
    io::OFastqPairedStream &unbinned_reads_ostream_;

    const Binning &binning_;
    const SoftBinsAssignment& soft_edge_labels_;
    const BinningAssignmentStrategy& assignment_strategy_;

    double bin_weight_threshold_;
};

void SplitAndWriteReads(const debruijn_graph::Graph &graph,
                        SequencingLib &lib,
                        const Binning &binning,
                        const SoftBinsAssignment& edge_soft_labels,
                        const BinningAssignmentStrategy& assignment_strategy,
                        const std::string &work_dir,
                        const std::string &prefix,
                        unsigned nthreads,
                        const double bin_weight_threshold) {
    if (!io::ReadConverter::LoadLibIfExists(lib)) {
        std::unique_ptr<ThreadPool::ThreadPool> pool;
        if (nthreads > 1)
            pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);
        io::ReadConverter::ConvertToBinary(lib, pool.get());
    }

    io::OFastqPairedStream unbinned_reads_ostream(fs::append_path(prefix, "unbinned_1.fastq"),
                                                  fs::append_path(prefix, "unbinned_2.fastq"),
                                                  lib.orientation());

    std::vector<io::OFastqPairedStream> binned_reads_ostreams;
    for (size_t bin_id = 0; bin_id < binning.bins().size(); ++bin_id) {
        auto bin_label = binning.bin_labels().at(bin_id);
        std::replace(bin_label.begin(), bin_label.end(), '/', '_');
        const std::pair<std::string, std::string> cur_part_paired_reads_filenames
                {fs::append_path(prefix, bin_label + "_1.fastq"),
                 fs::append_path(prefix, bin_label + "_2.fastq")};
        binned_reads_ostreams.emplace_back(cur_part_paired_reads_filenames.first,
                                           cur_part_paired_reads_filenames.second,
                                           lib.orientation());
    }

    SequenceMapperNotifier notifier;
    auto mapper = alignment::ShortKMerReadMapper(graph, work_dir);
    ReadsPathsToComponentsListener listener(binning,
                                            binned_reads_ostreams,
                                            unbinned_reads_ostream,
                                            edge_soft_labels,
                                            assignment_strategy,
                                            bin_weight_threshold);
    notifier.Subscribe(&listener);
    auto paired_streams = paired_binary_readers(lib, /*followed by rc*/ false, 0,
                                                /*include merged*/true);
    notifier.ProcessLibrary(paired_streams, mapper);
}

}


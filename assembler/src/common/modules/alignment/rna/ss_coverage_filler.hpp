//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <modules/alignment/sequence_mapper_notifier.hpp>
#include "modules/alignment/rna/ss_coverage.hpp"

namespace debruijn_graph {

class SSCoverageFiller: public SequenceMapperListener {
private:
    const Graph& g_;

    SSCoverageStorage& storage_;

    std::vector<SSCoverageStorage> tmp_storages_;

    bool symmetric_;

    void ProcessRange(size_t thread_index, const MappingPath<EdgeId>& read) {
        for (size_t i = 0; i < read.size(); ++i) {
            const auto& range = read[i].second;
            size_t kmer_count = range.mapped_range.end_pos - range.mapped_range.start_pos;
            tmp_storages_[thread_index].IncreaseKmerCount(read[i].first, kmer_count, symmetric_);
        }
    }

    DECL_LOGGER("SSCoverageFiller");
public:
    SSCoverageFiller(const Graph& g, SSCoverageStorage& storage, bool symmertic = false):
        g_(g), storage_(storage), tmp_storages_(), symmetric_(symmertic) {}

    void StartProcessLibrary(size_t threads_count) override {
        tmp_storages_.clear();

        for (size_t i = 0; i < threads_count; ++i) {
            tmp_storages_.emplace_back(g_);
        }
    }

    void StopProcessLibrary() override {
        for (auto& storage : tmp_storages_)
            storage.Clear();
        storage_.RecalculateCoverage();
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleRead& /* r */, const MappingPath<EdgeId>& read) override {
        ProcessRange(thread_index, read);
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleReadSeq& /* r */, const MappingPath<EdgeId>& read) override {
        ProcessRange(thread_index, read);
    }

    void MergeBuffer(size_t thread_index) override {
        for (const auto& it : tmp_storages_[thread_index])
            storage_.IncreaseKmerCount(it.first, size_t(it.second));
        tmp_storages_[thread_index].Clear();
    }
};


class SSBinCoverageFiller: public SequenceMapperListener {
private:
    SSCoverageSplitter& storage_;

    std::vector<SSCoverageSplitter> tmp_storages_;

    void ProcessRange(size_t thread_index, const MappingPath<EdgeId>& read) {
        for (size_t i = 0; i < read.size(); ++i) {
            auto range = read.mapping_at(i).mapped_range;
            tmp_storages_[thread_index].IncreaseKmerCount(read.edge_at(i), range);
        }
    }

    DECL_LOGGER("SSCoverageFiller");
public:
    explicit SSBinCoverageFiller(SSCoverageSplitter& storage):
        storage_(storage), tmp_storages_() {}

    void StartProcessLibrary(size_t threads_count) override {
        tmp_storages_.clear();

        for (size_t i = 0; i < threads_count; ++i) {
            tmp_storages_.emplace_back(storage_.g(), storage_.bin_size(), storage_.min_edge_len(),
                storage_.min_edge_coverage(), storage_.coverage_margin(), storage_.min_flanking_coverage());
        }
    }

    void StopProcessLibrary() override {
        for (auto& storage : tmp_storages_)
            storage.Clear();
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleRead& /* r */, const MappingPath<EdgeId>& read) override {
        ProcessRange(thread_index, read);
    }

    void ProcessSingleRead(size_t thread_index, const io::SingleReadSeq& /* r */, const MappingPath<EdgeId>& read) override {
        ProcessRange(thread_index, read);
    }

    void MergeBuffer(size_t thread_index) override {
        storage_.MergeOther(tmp_storages_[thread_index]);
        tmp_storages_[thread_index].Clear();
    }
};

}

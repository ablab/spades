//
// Created by andrey on 23.05.17.
//

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

}
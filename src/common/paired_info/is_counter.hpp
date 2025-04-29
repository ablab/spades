//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef IS_COUNTER_HPP_
#define IS_COUNTER_HPP_

#include "paired_info/insert_size_refiner.hpp"

#include "alignment/sequence_mapper_notifier.hpp"

namespace debruijn_graph {

using namespace omnigraph;

class InsertSizeCounter: public SequenceMapperListener {

public:

    InsertSizeCounter(const Graph &graph,
            size_t edge_length_threshold,
            bool ignore_negative = false)
        : graph_(graph),
          edge_length_threshold_(edge_length_threshold),
          ignore_negative_(ignore_negative) {
    }

    HistType hist() { return hist_; }
    size_t total() const { return total_.total_; }
    size_t mapped() const { return counted_.total_; }
    size_t negative() const { return negative_.total_; }


    void StartProcessLibrary(size_t threads_count) override {
        hist_.clear();
        tmp_hists_ = std::vector<HistType>(threads_count);

        total_ = count_data(threads_count);
        counted_ = count_data(threads_count);
        negative_ = count_data(threads_count);
    }

    void StopProcessLibrary() override {
        tmp_hists_.clear();
        total_.merge();
        counted_.merge();
        negative_.merge();
    }

    void ProcessPairedRead(size_t thread_index,
                           const io::PairedRead& r,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(thread_index, read1, read2, (int) r.second().size(),
                          (int) r.first().GetLeftOffset() + (int) r.second().GetRightOffset());
    }

    void ProcessPairedRead(size_t thread_index,
                           const io::PairedReadSeq& r,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(thread_index, read1, read2, (int) r.second().size(),
                          (int) r.first().GetLeftOffset() + (int) r.second().GetRightOffset());
    }

    void MergeBuffer(size_t thread_index) override {
        for (const auto& kv: tmp_hists_[thread_index])
            hist_[kv.first] += kv.second;

        tmp_hists_[thread_index].clear();
    }

    void FindMean(double& mean, double& delta, std::map<size_t, size_t>& percentiles) const {
        find_mean(hist_, mean, delta, percentiles);
    }

    void FindMedian(double& median, double& mad, HistType& histogram) const {
        find_median(hist_, median, mad, histogram);
    }

private:
    void ProcessPairedRead(size_t thread_index,
                           const MappingPath<EdgeId>& read1, const MappingPath<EdgeId>& read2,
                           int read2_size,
                           int is_delta) {

        ++total_.arr_[thread_index];

        if (read1.size() == 1 && read2.size() == 1 &&
            read2.simple_path().front() == read1.simple_path().front() &&
            graph_.length(read1.simple_path().front()) >= edge_length_threshold_) {

            auto mapping_edge_1 = read1.front().second;
            auto mapping_edge_2 = read2.front().second;

            int read1_start = (int) mapping_edge_1.mapped_range.start_pos - (int) mapping_edge_1.initial_range.start_pos ;
            TRACE("Read 1: " << (int) mapping_edge_1.mapped_range.start_pos << " - " << (int) mapping_edge_1.initial_range.start_pos << " = " << read1_start);
            int read2_start = (int) mapping_edge_2.mapped_range.start_pos - (int) mapping_edge_2.initial_range.start_pos;
            TRACE("Read 2: " << (int) mapping_edge_2.mapped_range.start_pos << " - " << (int) mapping_edge_2.initial_range.start_pos << " = " << read2_start);
            int is = read2_start - read1_start + read2_size + is_delta;
            TRACE("IS: " << read2_start << " - " <<  read1_start << " + " << (int) is_delta << " = " << is);

            if (is > 0 || ignore_negative_) {
                tmp_hists_[thread_index][is] += 1;
                ++counted_.arr_[thread_index];
            } else {
                ++negative_.arr_[thread_index];
            }

        }

    }
    struct count_data {
        size_t total_;
        std::vector<size_t> arr_;
        count_data()
            : total_(0) {}

        count_data(size_t nthreads)
            : total_(0), arr_(nthreads, 0) {}

        void inc(size_t i) { ++arr_[i]; }
        void merge() {
            for (size_t i = 0; i < arr_.size(); ++i) {
                total_ += arr_[i];
                arr_[i] = 0;
            }
        }

        void BinWrite(std::ostream &os) const {
            io::binary::BinWrite(os, total_);
        }

        void BinRead(std::istream &is) {
            io::binary::BinRead(is, total_);
        }
    };

    void Serialize(std::ostream &os) const override {
        total_.BinWrite(os);
        counted_.BinWrite(os);
        negative_.BinWrite(os);
        io::binary::BinWrite(os, hist_);
    }

    void Deserialize(std::istream &is) override {
        total_.BinRead(is);
        counted_.BinRead(is);
        negative_.BinRead(is);
        io::binary::BinRead(is, hist_);
    }

    void MergeFromStream(std::istream &is) override {
        InsertSizeCounter remote(*this);
        remote.Deserialize(is);
        total_.total_ += remote.total_.total_;
        counted_.total_ += remote.counted_.total_;
        negative_.total_ += remote.negative_.total_;

        for (const auto& kv: remote.hist_) {
            hist_[kv.first] += kv.second;
        }
    }

private:
    const Graph &graph_;

    HistType hist_;
    std::vector<HistType> tmp_hists_;

    count_data total_;
    count_data counted_;
    count_data negative_;

    size_t edge_length_threshold_;
    bool ignore_negative_;
};

}


#endif /* IS_COUNTER_HPP_ */

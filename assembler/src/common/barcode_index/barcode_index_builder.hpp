//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index.hpp"

#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/read_processor.hpp"
#include "alignment/sequence_mapper_notifier.hpp"
#include "alignment/sequence_mapper.hpp"

#include <mutex>
#include <random>
#include <shared_mutex>
#include <unordered_set>

namespace barcode_index {

//todo templatize
class ConcurrentBufferFiller {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;
    using FrameConcurrentBuffer = FrameConcurrentBarcodeIndexBuffer<Graph>;
    using SequenceMapper = debruijn_graph::SequenceMapper<Graph>;
    using MappingPath = omnigraph::MappingPath<EdgeId>;
    using MappingRange = omnigraph::MappingRange;

    class BarcodeEncoder {
      public:
        BarcodeEncoder(size_t start):
            codes_(),
            start_(start)
        { }

        BarcodeId add(const std::string &barcode) {
            std::unique_lock<std::shared_timed_mutex> lock(mutex_);
            size_t encoder_size = codes_.size();
            codes_[barcode] = encoder_size + start_;
            return encoder_size;
        }

        auto begin() const {
            return codes_.begin();
        }
        auto end() {
            std::shared_lock<std::shared_timed_mutex> lock(mutex_);
            return codes_.end();
        }
        auto find(const string& barcode) {
            std::shared_lock<std::shared_timed_mutex> lock(mutex_);
            return codes_.find(barcode);
        }
        size_t size() const {
            return codes_.size();
        }
      private:
        std::unordered_map <std::string, BarcodeId> codes_;
        std::shared_timed_mutex mutex_;
        const size_t start_;
    };

    ConcurrentBufferFiller(const Graph &g,
                           FrameConcurrentBuffer &buf,
                           const SequenceMapper &mapper,
                           const std::vector<std::string> &barcode_prefices,
                           size_t starting_barcode, bool is_tellseq) : g_(g),
                                                                       buf_(buf),
                                                                       mapper_(mapper),
                                                                       barcode_prefices_(barcode_prefices),
                                                                       encoder_(starting_barcode),
                                                                       is_tellseq_(is_tellseq) { Init(); }

    bool operator()(std::unique_ptr<io::PairedRead> r) {
        const Sequence& read1 = r->first().sequence();
        const Sequence& read2 = r->second().sequence();
        MappingPath path1 = mapper_.MapSequence(read1);
        MappingPath path2 = mapper_.MapSequence(read2);
        auto barcode = GetBarcode(std::move(r));
        if (barcode.empty()) {
            TRACE("Empty barcode")
            return false;
        }
        if (path1.size() > 0 and path2.size() > 0) {
            TRACE("non-empty pair");
        }

        ProcessPairedRead(barcode, path1, path2);
        return false;
    }

    size_t GetNumberOfBarcodes() const {
        return encoder_.size();
    }

    FrameConcurrentBuffer& GetBuffer() {
        return buf_;
    }

  private:
    void ProcessPairedRead(const std::string &barcode_string,
                           const MappingPath& path1,
                           const MappingPath& path2) {

        auto code_result = encoder_.find(barcode_string);
        BarcodeId barcode;
        if (code_result == encoder_.end()) {
            barcode = encoder_.add(barcode_string);
        } else {
            barcode = code_result->second;
        }
        InsertMappingPath(barcode, path1);
        InsertMappingPath(barcode, path2);
    }
    void Init() {
        buf_.InitialFillMap();
    }

    std::string GetBarcode(std::unique_ptr<io::PairedRead> r) {
        if (not is_tellseq_) {
            auto left_barcode_string = GetTenXBarcodeFromRead(r->first().name(), barcode_prefices_);
            auto right_barcode_string = GetTenXBarcodeFromRead(r->second().name(), barcode_prefices_);

            if (left_barcode_string.empty() or left_barcode_string != right_barcode_string) {
                TRACE(left_barcode_string);
                TRACE(right_barcode_string);
                std::string empty;
                return empty;
            }
            return left_barcode_string;
        } else {
            io::PairedRead* paired = r.get();
            return static_cast<io::TellSeqRead*>(paired)->aux().sequence().str();
        }
    }

    string GetTenXBarcodeFromRead(const std::string &read_name, const std::vector<string>& barcode_prefixes) {
        for (const auto& prefix: barcode_prefixes) {
            size_t prefix_len = prefix.size();
            size_t start_pos = read_name.find(prefix);
            if (start_pos != string::npos) {
                string barcode = GetBarcodeFromStartPos(start_pos + prefix_len, read_name);
                TRACE(barcode);
                return barcode;
            }
        }
        return "";
    }
    string GetBarcodeFromStartPos(const size_t start_pos, const string& read_id) {
        string result = "";
        for (auto it = read_id.begin() + start_pos; it != read_id.end(); ++it) {
            if (not is_nucl(*it)) {
                return result;
            }
            result.push_back(*it);
        }
        return result;
    }

    void InsertMappingPath(BarcodeId &barcode, const MappingPath &path) {
        for (size_t i = 0; i < path.size(); i++) {
            //todo restore tail threshold if needed
            EdgeId edge = path[i].first;
            const auto &range = path[i].second.mapped_range;
            buf_.InsertBarcode(barcode, edge, 1, range);
            buf_.InsertBarcode(barcode, g_.conjugate(edge), 1, range.Invert(g_.length(edge)));
        }
    }

    const Graph& g_;
    FrameConcurrentBuffer &buf_;
    const SequenceMapper &mapper_;
    const std::vector<std::string> barcode_prefices_;
    BarcodeEncoder encoder_;
    bool is_tellseq_;

    DECL_LOGGER("ConcurrentBufferFiller");
};

class FrameBarcodeIndexBuilder {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;
    using SequenceMapper = debruijn_graph::SequenceMapper<Graph>;

    FrameBarcodeIndexBuilder(const Graph &g,
                             const SequenceMapper &mapper,
                             const std::vector<std::string> &barcode_prefices,
                             size_t frame_size,
                             size_t num_threads) :
        g_(g),
        mapper_(mapper),
        barcode_prefices_(barcode_prefices),
        frame_size_(frame_size),
        num_threads_(num_threads) {}

    template<class ReadType>
    void ConstructBarcodeIndex(io::ReadStreamList<ReadType> read_streams,
                               FrameBarcodeIndex<Graph> &barcode_index,
                               const io::SequencingLibraryBase &lib,
                               bool is_tellseq);

    void DownsampleBarcodeIndex(FrameBarcodeIndex<Graph> &downsampled_index, FrameBarcodeIndex<Graph> &original_index, double sampling_factor) {
        std::unordered_set<BarcodeId> barcodes;
        for (auto it = original_index.begin(); it != original_index.end(); ++it) {
            const auto &barcode_distribution = it->second.GetDistribution();
            for (const auto &entry: barcode_distribution) {
                BarcodeId current_barcode = entry.first;
                barcodes.insert(current_barcode);
            }
        }
        INFO("Number of encountered barcodes: " << barcodes.size());
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distr(.0, 1.0);
        std::unordered_set<BarcodeId> passed_barcodes;
        for (const auto &barcode: barcodes) {
            if (math::le(distr(gen), sampling_factor)) {
                passed_barcodes.insert(barcode);
            }
        }
        INFO("Passed barcodes: " << passed_barcodes.size());

        downsampled_index.InitialFillMap();
        auto barcode_filter = [&passed_barcodes](const auto &barcode_entry) {
          return passed_barcodes.find(barcode_entry.first) != passed_barcodes.end();
        };
        size_t final_barcodes = 0;
        for (auto it = original_index.begin(); it != original_index.end(); ++it) {
            auto &to = downsampled_index.edge_to_entry_[it->first].barcode_distribution_;
            auto &from = it->second.barcode_distribution_;
            std::copy_if(std::make_move_iterator(from.begin()), std::make_move_iterator(from.end()),
                         std::inserter(to, to.end()), barcode_filter);
        }
    }

  private:
    const Graph& g_;
    const SequenceMapper &mapper_;
    const std::vector<std::string> barcode_prefices_;
    size_t frame_size_;
    size_t num_threads_;
};
template<class ReadType>
void FrameBarcodeIndexBuilder::ConstructBarcodeIndex(io::ReadStreamList<ReadType> read_streams,
                                                     FrameBarcodeIndex<Graph> &barcode_index,
                                                     const io::SequencingLibraryBase &lib,
                                                     bool is_tellseq) {
    {
        size_t starting_barcode = 0;
        size_t counter = 0;
        barcode_index.SetFrameSize(frame_size_);
        barcode_index.InitialFillMap();
        for (auto &stream: read_streams) {
            DEBUG("Processing stream " << counter << " , currently " << starting_barcode << " barcodes");
            FrameConcurrentBarcodeIndexBuffer<debruijn_graph::Graph> buffer(g_, frame_size_);
            ConcurrentBufferFiller buffer_filler(g_, buffer, mapper_, barcode_prefices_, starting_barcode, is_tellseq);
            hammer::ReadProcessor read_processor(static_cast<unsigned int>(num_threads_));
            read_processor.Run(stream, buffer_filler);
            starting_barcode += buffer_filler.GetNumberOfBarcodes();
            DEBUG("Update");
            barcode_index.Update(buffer);
            DEBUG("Finished update");
        }
        INFO(starting_barcode << " total barcodes in the barcode index");
    }

}
} //namespace barcode index

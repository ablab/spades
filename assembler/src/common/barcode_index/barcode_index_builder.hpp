//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index.hpp"

#include "modules/alignment/sequence_mapper_notifier.hpp"

#include <mutex>
#include <shared_mutex>

namespace barcode_index {

//todo templatize
class FrameMapperBuilder: public debruijn_graph::SequenceMapperListener {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;
    using FrameBarcodeMapper = FrameBarcodeIndex<Graph>;
    typedef omnigraph::MappingPath<EdgeId> MappingPath;
    typedef omnigraph::MappingRange MappingRange;

    class BarcodeEncoder {
      public:
        BarcodeEncoder():
            codes_()
        { }

        BarcodeId add(const std::string &barcode) {
            std::unique_lock<std::shared_timed_mutex> lock(mutex_);
            size_t encoder_size = codes_.size();
            codes_[barcode] = encoder_size;
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
    };

    FrameMapperBuilder(const Graph &g,
                       FrameBarcodeMapper &index,
                       size_t num_threads,
                       size_t frame_size,
                       const std::vector<std::string> &barcode_prefices) : g_(g),
                                                                           buf_(num_threads,
                                                                                FrameBarcodeMapper(g, frame_size)),
                                                                           data_(index),
                                                                           frame_size_(frame_size),
                                                                           barcode_prefices_(barcode_prefices),
                                                                           encoder_() {Init();}

    void MergeBuffer(size_t i) override {
        for (auto it = buf_[i].begin(); it != buf_[i].end(); ++it) {
            EdgeId edge = it->first;
            const auto &entry = it->second;
            //todo might be faster to iterate over the main storage if a buffer is large enough
            for (auto barcode_it = entry.begin(); barcode_it != entry.end(); ++barcode_it) {
                data_.edge_to_entry_[edge].InsertInfo(barcode_it->first, barcode_it->second);
            }
        }
        DEBUG("Merged buffer");
    }

    void ProcessPairedRead(size_t idx,
                           const io::PairedRead& r,
                           const MappingPath& read1,
                           const MappingPath& read2) override {
        ProcessPairedRead(buf_[idx], r.first().name(), r.second().name(), read1, read2);
    }

    size_t GetNumberOfBarcodes() const {
        return encoder_.size();
    }

//    void ExportToIndex(HicCovPlotIndex& index);

  private:
    void ProcessPairedRead(FrameBarcodeMapper &buf,
                           const std::string &name1,
                           const std::string &name2,
                           const MappingPath& path1,
                           const MappingPath& path2) {
        std::string barcode_string = GetTenXBarcodeFromRead(name1, barcode_prefices_);
        std::string second_barcode_string = GetTenXBarcodeFromRead(name2, barcode_prefices_);
        if (barcode_string.empty() or barcode_string != second_barcode_string) {
            return;
        }

        auto code_result = encoder_.find(barcode_string);
        BarcodeId barcode;
        if (code_result == encoder_.end()) {
            barcode = encoder_.add(barcode_string);
        } else {
            barcode = code_result->second;
        }
        InsertMappingPath(buf, barcode, path1);
        InsertMappingPath(buf, barcode, path2);
    }
    void Init() {
        data_.InitialFillMap();
        for (auto &buf: buf_) {
            buf.InitialFillMap();
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
            if (std::isspace(*it)) {
                return result;
            }
            result.push_back(*it);
        }
        return result;
    }

    void InsertMappingPath(FrameBarcodeMapper &buf, BarcodeId &barcode, const MappingPath &path) {
        for (size_t i = 0; i < path.size(); i++) {
            //todo restore tail threshold if needed
            EdgeId edge = path[i].first;
            const auto &range = path[i].second.mapped_range;
            InsertBarcode(buf, barcode, edge, 1, range);
            InsertBarcode(buf, barcode, g_.conjugate(edge), 1, range.Invert(g_.length(edge)));
        }
    }

    void InsertBarcode(FrameBarcodeMapper &buf, const BarcodeId &barcode, const EdgeId &edge, size_t count, const Range &range) {
        VERIFY_DEV(buf.edge_to_entry_.find(edge) != buf.edge_to_entry_.end());
        buf.edge_to_entry_.at(edge).InsertBarcode(barcode, count, range);
    }

    const Graph& g_;
    std::vector<FrameBarcodeMapper> buf_;
    FrameBarcodeMapper &data_;
    size_t frame_size_;
    const std::vector<std::string> barcode_prefices_;
    BarcodeEncoder encoder_;

    DECL_LOGGER("FrameMapperBuilder");
};
}
//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index.hpp"

#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/read_processor.hpp"
#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include <mutex>
#include <shared_mutex>

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

    ConcurrentBufferFiller(const Graph &g,
                           FrameConcurrentBuffer &buf,
                           const SequenceMapper &mapper,
                           const std::vector<std::string> &barcode_prefices) : g_(g),
                                                                               buf_(buf),
                                                                               mapper_(mapper),
                                                                               barcode_prefices_(barcode_prefices),
                                                                               encoder_() { Init(); }

    bool operator()(std::unique_ptr<io::PairedRead> r) {
        const Sequence& read1 = r->first().sequence();
        const Sequence& read2 = r->second().sequence();
        MappingPath path1 = mapper_.MapSequence(read1);
        MappingPath path2 = mapper_.MapSequence(read2);

        ProcessPairedRead(r->first().name(), r->second().name(), path1, path2);
        return false;
    }

    size_t GetNumberOfBarcodes() const {
        return encoder_.size();
    }

    FrameConcurrentBuffer& GetBuffer() {
        return buf_;
    }

  private:
    void ProcessPairedRead(const std::string &name1,
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
        InsertMappingPath(barcode, path1);
        InsertMappingPath(barcode, path2);
    }
    void Init() {
        buf_.InitialFillMap();
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

    DECL_LOGGER("ConcurrentBufferFiller");
};

class FrameBarcodeIndexBuilder {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;

    FrameBarcodeIndexBuilder(ConcurrentBufferFiller &buffer_filler, size_t num_threads) : buffer_filler_(buffer_filler),
                                                                                          num_threads_(num_threads) {}

    void ConstructBarcodeIndex(FrameBarcodeIndex<Graph> &barcode_index, const io::SequencingLibraryBase &lib) {
        auto read_streams = io::paired_easy_readers(lib, false, 0);
        for (auto &stream: read_streams) {
            hammer::ReadProcessor read_processor(static_cast<unsigned int>(num_threads_));
            read_processor.Run(stream, buffer_filler_);
            auto &buffer = buffer_filler_.GetBuffer();
            barcode_index.SetFrameSize(buffer.GetFrameSize());
            barcode_index.MoveAssign(buffer);
        }
    }
  private:
    ConcurrentBufferFiller &buffer_filler_;
    size_t num_threads_;
};
} //namespace barcode index
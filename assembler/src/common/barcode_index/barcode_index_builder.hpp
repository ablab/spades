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
                           size_t starting_barcode) : g_(g),
                                                      buf_(buf),
                                                      mapper_(mapper),
                                                      barcode_prefices_(barcode_prefices),
                                                      encoder_(starting_barcode) { Init(); }

    bool operator()(std::unique_ptr<io::PairedRead> r) {
        const Sequence& read1 = r->first().sequence();
        const Sequence& read2 = r->second().sequence();
        MappingPath path1 = mapper_.MapSequence(read1);
        MappingPath path2 = mapper_.MapSequence(read2);
        if (path1.size() > 0 and path2.size() > 0) {
            TRACE("non-empty pair");
        }

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
            TRACE(barcode_string);
            TRACE(second_barcode_string);
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

    void ConstructBarcodeIndex(FrameBarcodeIndex<Graph> &barcode_index, const io::SequencingLibraryBase &lib) {
        auto read_streams = io::paired_easy_readers(lib, false, 0);
        size_t starting_barcode = 0;
        size_t counter = 0;
        barcode_index.SetFrameSize(frame_size_);
        barcode_index.InitialFillMap();
        for (auto &stream: read_streams) {
            INFO("Processing stream " << counter << " , currently " << starting_barcode << " barcodes");
            FrameConcurrentBarcodeIndexBuffer<debruijn_graph::Graph> buffer(g_, frame_size_);
            ConcurrentBufferFiller buffer_filler(g_, buffer, mapper_, barcode_prefices_, starting_barcode);
            hammer::ReadProcessor read_processor(static_cast<unsigned int>(num_threads_));
            read_processor.Run(stream, buffer_filler);
            starting_barcode += buffer_filler.GetNumberOfBarcodes();
            INFO("Update");
            barcode_index.Update(buffer);
            INFO("Finished update");
        }
    }
  private:
    const Graph& g_;
    const SequenceMapper &mapper_;
    const std::vector<std::string> barcode_prefices_;
    size_t frame_size_;
    size_t num_threads_;
};
} //namespace barcode index
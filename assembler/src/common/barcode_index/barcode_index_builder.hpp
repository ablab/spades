//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/edge_index.hpp"
#include "modules/alignment/kmer_mapper.hpp"
#include "modules/alignment/sequence_mapper.hpp"

namespace barcode_index {

    template<class Graph, class BarcodeEntryT>
    class BarcodeIndexBuilder {
    public:
        typedef debruijn_graph::EdgeIndex<Graph> Index;
        typedef debruijn_graph::KmerMapper<Graph> KmerSubs;
        typedef debruijn_graph::KmerFreeEdgeIndex<Graph, utils::DefaultStoring> InnerIndex;
        typedef typename InnerIndex::KeyWithHash KeyWithHash;
        typedef typename barcode_index::BarcodeIndex<Graph, BarcodeEntryT> BarcodeIndexT;
        typedef typename Graph::EdgeId EdgeId;

        BarcodeIndexBuilder(BarcodeIndexT &index, size_t tail_threshold) :
                g_(index.GetGraph()),
                index_(index),
                tail_threshold_(tail_threshold),
                barcode_codes_() {}

        ~BarcodeIndexBuilder() {}

        BarcodeIndexT GetMapper() {
            return index_;
        }

        void FillMapFromDemultiplexedDataset(const Index &index, const KmerSubs &kmer_mapper,
                                             const std::string &tslr_dataset) {
            auto lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph,Index>>(g_, index, kmer_mapper);
            //Process every barcode from truspades dataset
            for (size_t i = 0; i < lib_vec.size(); ++i) {
                std::string barcode_string = lib_vec[i].barcode_;
                uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                BarcodeId barcode(barcode_int);

                std::shared_ptr<io::ReadStream<io::PairedRead>> paired_stream =
                        std::make_shared<io::SeparatePairedReadStream>(lib_vec[i].left_, lib_vec[i].right_, 1);
                io::PairedRead read;
                while (!paired_stream->eof()) {
                    *paired_stream >> read;
                    auto path_first = mapper->MapRead(read.first());
                    auto path_second = mapper->MapRead(read.second());
                    InsertMappingPath(barcode, path_first);
                    InsertMappingPath(barcode, path_second);
                }
                VERBOSE_POWER_T2(i, 100,
                                 "Processed " << i << " barcodes from " << lib_vec.size() << " ("
                                              << i * 100 / lib_vec.size()
                                              << "%)");
                if (lib_vec.size() > 10 && i % (lib_vec.size() / 10 + 1) == 0) {
                    INFO("Processed " << i << " barcodes from " << lib_vec.size() << " (" << i * 100 / lib_vec.size()
                                      << "%)");
                }
            }
        }

        void FillMapUsingKmerMultisetParallel(const Index &index,
                                              const KmerSubs &kmer_mapper,
                                              const std::string &tslr_dataset,
                                              size_t n_threads) {
            const auto &lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index>>(g_, index, kmer_mapper);
            const auto &bucket_vec = SplitLibrary(lib_vec, n_threads);
            INFO("Library splitted");
            for (size_t i = 0; i < bucket_vec.size(); ++i) {
#pragma omp parallel for num_threads(n_threads)
                for (size_t j = 0; j < bucket_vec[i].size(); ++j) {
                    const auto &lib = bucket_vec[i][j];
                    std::string barcode_string = lib.barcode_;
                    uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                    BarcodeId barcode(barcode_int);
                    ReadStreamPtr paired_stream = std::make_shared<io::SeparatePairedReadStream>(lib.left_,
                                                                                                 lib.right_, 1);
                    const KmerMultiset &kmer_multiset = BuildKmerMultisetFromStream(paired_stream);
                    for (auto it = kmer_multiset.cbegin(); it != kmer_multiset.cend(); ++it) {
                        size_t count = it->second;
                        const auto &edge_and_offset = index.get(it->first);
                        EdgeId edge = edge_and_offset.first;
                        size_t offset = edge_and_offset.second;
                        if (edge.int_id() != 0) {
#pragma omp critical
                            {
                                InsertBarcodeWithRange(barcode, edge, Range(offset, offset + 1), count);
                            }
                        }
                    }
                }
                INFO((i + 1) * n_threads << " barcodes processed.");
            }
        }

        template <class Mapper>
        void FillMapFrom10XReads(std::vector<io::SingleStream> &reads, std::shared_ptr<Mapper> read_mapper) {
            INFO("Starting barcode index construction from 10X reads")
            //Process every read from 10X dataset
            io::SingleRead read;
            size_t counter = 0;
            const std::vector<string> barcode_prefixes = {"BC:Z:", "BX:Z:"};
            for (auto &stream: reads) {
                while (!stream.eof()) {
                    stream >> read;
                    string barcode_string = GetTenXBarcodeFromRead(read, barcode_prefixes);

                    if (barcode_string != "") {
                        barcode_codes_.AddBarcode(barcode_string);
                        uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                        BarcodeId barcode(barcode_int);
                        const auto &path = read_mapper->MapRead(read);
                        InsertMappingPath(barcode, path);
                    }
                    counter++;
                    VERBOSE_POWER_T2(counter, 100, "Processed " << counter << " reads.");
                }
            }
            index_.SetNumberOfBarcodes(barcode_codes_.GetSize());
            INFO("FillMap finished")
        }

        void FillMap(std::vector<io::SingleStream> &reads) {
            InitialFillMap();
            auto mapper = std::make_shared<alignment::BWAReadMapper<Graph>>(g_);
            FillMapFrom10XReads(reads, mapper);
            return;
        }

        void FillMap(std::vector<io::SingleStream> &reads, const Index &index, const KmerSubs &kmer_mapper) {
            InitialFillMap();
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index>>(g_, index, kmer_mapper);
            FillMapFrom10XReads(reads, mapper);
            return;
        }
        virtual void InitialFillMap() = 0;
        DECL_LOGGER("BarcodeMapperBuilder");

    protected:
        typedef std::shared_ptr<io::ReadStream<io::PairedRead>> ReadStreamPtr;

        void InsertEntry(const EdgeId &edge, BarcodeEntryT &entry) {
            index_.InsertEntry(edge, entry);
        }

        string GetTenXBarcodeFromRead(const io::SingleRead &read, const std::vector<string>& barcode_prefixes) {
            for (const auto& prefix: barcode_prefixes) {
                size_t prefix_len = prefix.size();
                size_t start_pos = read.name().find(prefix);
                if (start_pos != string::npos) {
                    string barcode = GetBarcodeFromStartPos(start_pos + prefix_len, read.name());
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

        KmerMultiset BuildKmerMultisetFromStream(ReadStreamPtr read_stream) {
            KmerMultiset kmer_multiset;
            size_t read_counter = 0;
            io::PairedRead read;
            size_t kmer_len = g_.k() + 1;
            while (!read_stream->eof()) {
                read_counter += 2;
                *read_stream >> read;
                ExtractKmersFromSeq(kmer_multiset, read.first(), kmer_len);
                ExtractKmersFromSeq(kmer_multiset, read.second(), kmer_len);
            }
            return kmer_multiset;
        }

        void ExtractKmersFromSeq(KmerMultiset &kmer_multiset, const io::SingleRead &read, size_t kmer_len) {
            if (read.IsValid()) {
                const Sequence &seq = read.sequence();
                Kmer kmer = seq.start<Kmer>(kmer_len);
                for (size_t j = kmer_len; j < seq.size(); ++j) {
                    kmer_multiset.Insert(kmer);
                    kmer <<= seq[j];
                }
                kmer_multiset.Insert(kmer);
            }
        }


        void InsertBarcode(const BarcodeId &barcode, const EdgeId &edge, size_t count, const Range &range) {
            VERIFY_DEV(index_.edge_to_entry_.find(edge) != index_.edge_to_entry_.end());
            index_.edge_to_entry_.at(edge).InsertBarcode(barcode, count, range);
        }

        bool IsAtEdgeTail(const EdgeId &edge, const Range &range) {
            return range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool IsAtEdgeHead(const Range &range) {
            return range.end_pos < tail_threshold_;
        }

        void InsertMappingPath(const BarcodeId &barcode, const MappingPath <EdgeId> &path) {
            for (size_t j = 0; j < path.size(); j++) {
                InsertBarcodeWithRange(barcode, path[j].first, path[j].second.mapped_range, 1);
            }
        }

        void InsertBarcodeWithRange(const BarcodeId &barcode, const EdgeId &edge,
                                    const Range &range, size_t count) {
            if (IsAtEdgeHead(range))
                InsertBarcode(barcode, edge, count, range);
            if (IsAtEdgeTail(edge, range))
                InsertBarcode(barcode, g_.conjugate(edge), count, range.Invert(g_.length(edge)));
        }


        std::vector <tslr_barcode_library> GetLibrary(const string &tslr_dataset_name) {
            std::vector <tslr_barcode_library> lib_vec;
            std::ifstream fin;
            fin.open(tslr_dataset_name);
            string line;
            while (getline(fin, line)) {
                if (!line.empty()) {
                    istringstream tmp_stream(line);
                    tslr_barcode_library lib;
                    tmp_stream >> lib.barcode_;
                    tmp_stream >> lib.left_;
                    tmp_stream >> lib.right_;
                    barcode_codes_.AddBarcode(lib.barcode_);
                    lib_vec.push_back(lib);
                }
            }
            return lib_vec;
        }

        std::vector <std::vector<tslr_barcode_library>> SplitLibrary(const std::vector <tslr_barcode_library> &lib_vec,
                                                                     size_t bucket_size) {
            size_t lib_size = lib_vec.size();
            size_t buckets = lib_size / bucket_size;
            std::vector <std::vector<tslr_barcode_library>> bucket_vec(buckets);
            for (size_t i = 0; i < buckets; ++i) {
                for (size_t j = 0; j < bucket_size; ++j) {
                    bucket_vec[i].push_back(lib_vec[i * bucket_size + j]);
                }
            }
            std::vector <tslr_barcode_library> last_bucket;
            size_t last_elem = (lib_size / bucket_size) * bucket_size;
            for (size_t i = 0; i < lib_size % bucket_size; ++i) {
                last_bucket.push_back(lib_vec[last_elem + i]);
            }
            if (last_bucket.size() > 0) {
                bucket_vec.push_back(last_bucket);
            }
            size_t sum_size = 0;
            for (const auto &vec : bucket_vec) {
                sum_size += vec.size();
            }
            VERIFY(lib_size == sum_size);
            return bucket_vec;
        }

        const Graph &g_;
        BarcodeIndexT &index_;
        size_t tail_threshold_;
        BarcodeEncoder barcode_codes_;
    };

    template<class Graph>
    class FrameMapperBuilder : public BarcodeIndexBuilder<Graph, FrameEdgeEntry<Graph>> {
    public:
        typedef typename Graph::EdgeId EdgeId;
        typedef typename barcode_index::FrameBarcodeIndex<Graph> FrameBarcodeIndexT;

        FrameMapperBuilder(FrameBarcodeIndexT &index, const size_t tail_threshold, const size_t frame_size) :
                BarcodeIndexBuilder<Graph, FrameEdgeEntry<Graph>>(index, tail_threshold),
                frame_size_(frame_size) {}

        void InitialFillMap() override {
            FrameBarcodeIndexT& frame_index = dynamic_cast<FrameBarcodeIndexT&>(index_);
            frame_index.SetFrameSize(frame_size_);
            frame_index.InitialFillMap();
        }
     protected:
        using BarcodeIndexBuilder<Graph, FrameEdgeEntry<Graph>>::g_;
        using BarcodeIndexBuilder<Graph, FrameEdgeEntry<Graph>>::index_;
        const size_t frame_size_;
    };
}
#pragma once

#include "barcode_index.hpp"
#include "common/pipeline/library.hpp"

namespace barcode_index {
    template<class barcode_entry_t>
    class BarcodeIndexBuilder {
    protected:
        const Graph &g_;
        shared_ptr <BarcodeIndex<barcode_entry_t>> mapper_;
        size_t tail_threshold_;
        BarcodeEncoder barcode_codes_;
        typedef vector<io::SequencingLibrary<debruijn_graph::config::DataSetData>> lib_vector_t;

    public:
        BarcodeIndexBuilder(const Graph &g, size_t tail_threshold) :
                g_(g),
                mapper_(make_shared<BarcodeIndex<barcode_entry_t>>(g)),
                tail_threshold_(tail_threshold),
                barcode_codes_() {}

        ~BarcodeIndexBuilder() {}

        DECL_LOGGER("BarcodeMapperBuilder")

        shared_ptr <BarcodeIndex<barcode_entry_t>> GetMapper() {
            return mapper_;
        }

        void FillMapFromDemultiplexedDataset(const Index &index, const KmerSubs &kmer_mapper) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            auto lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared < debruijn_graph::BasicSequenceMapper < Graph, Index> >
            (g_, index, kmer_mapper);

            //Process every barcode from truspades dataset
            for (size_t i = 0; i < lib_vec.size(); ++i) {
                std::string barcode_string = lib_vec[i].barcode_;
                uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                BarcodeId barcode(barcode_int);

                std::shared_ptr <io::ReadStream<io::PairedRead>> paired_stream =
                        make_shared<io::SeparatePairedReadStream>(lib_vec[i].left_, lib_vec[i].right_, 1);
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
                ++mapper_->barcodes_number_;
            }
        }

        void FillMapUsingKmerMultisetParallel(const Index &index, const KmerSubs &kmer_mapper, size_t n_threads) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            const auto &lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared < debruijn_graph::BasicSequenceMapper < Graph, Index> >
            (g_, index, kmer_mapper);
            const auto &bucket_vec = SplitLibrary(lib_vec, n_threads);
            INFO("Library splitted");
            for (size_t i = 0; i < bucket_vec.size(); ++i) {
#pragma omp parallel for num_threads(n_threads)
                for (size_t j = 0; j < bucket_vec[i].size(); ++j) {
                    const auto &lib = bucket_vec[i][j];
                    std::string barcode_string = lib.barcode_;
                    uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                    BarcodeId barcode(barcode_int);

                    std::shared_ptr <io::ReadStream<io::PairedRead>> paired_stream =
                            make_shared<io::SeparatePairedReadStream>(lib.left_, lib.right_, 1);
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

        void FillMapUsingSubIndex(const Index &index, const KmerSubs &kmer_mapper) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            auto lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared < debruijn_graph::BasicSequenceMapper < Graph, Index> >
            (g_, index, kmer_mapper);
            size_t global_counter = 0;
            size_t counter = 0;

            //Process every barcode from truspades dataset
            for (size_t i = 0; i < lib_vec.size(); ++i) {
                std::string barcode_string = lib_vec[i].barcode_;
                uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                BarcodeId barcode(barcode_int);

                Index barcode_subindex(g_, cfg::get().tmp_dir);
                InnerIndex subindex = InnerIndex(g_, cfg::get().tmp_dir);

                io::ReadStreamList <io::SingleRead> streams;
                streams.push_back(io::EasyStream(lib_vec[i].left_, false /*followed_by_rc*/));
                streams.push_back(io::EasyStream(lib_vec[i].right_, false /*followed_by_rc*/));
                IndexBuilder().BuildIndexFromStream(subindex, streams, 0);

                INFO("Extracting kmer multiset using subindex")

                for (auto it = subindex.kmer_begin(); it.good(); ++it) {
                    Kmer kmer = Kmer(g_.k() + 1, *it);
                    KeyWithHash kh = subindex.ConstructKWH(kmer);
                    const auto &edge_and_offset = index.get(kmer);
                    EdgeId edge = edge_and_offset.first;
                    size_t offset = edge_and_offset.second;
                    size_t count = 0;
                    if (subindex.valid(kh)) {
                        count = subindex.get_value(kh).count;
                    }
                    if (edge.int_id() != 0 and count > 0) {
                        InsertBarcodeWithRange(barcode, edge, Range(offset, offset + 1), count);
                        ++counter;
                    }
                    ++global_counter;
                }
                ++mapper_->barcodes_number_;
            }
        }

        void FillMapFrom10XReads(const lib_vector_t& libs_10x) {
            INFO("Starting barcode index construction from 10X reads")
            auto mapper = std::make_shared < alignment::BWAReadMapper < Graph > >
                          (g_);

            auto streams = GetStreamsFromLibs(libs_10x);
            //Process every read from 10X dataset
            io::SingleRead read;
            size_t counter = 0;
            for (auto stream: streams) {
                while (!stream->eof()) {
                    *stream >> read;
                    string barcode_string = GetTenXBarcodeFromRead(read);

                    if (barcode_string != "") {
                        barcode_codes_.AddBarcode(barcode_string);
                        uint64_t barcode_int = barcode_codes_.GetCode(barcode_string);
                        BarcodeId barcode(barcode_int);

                        const auto &path = mapper->MapRead(read);
                        InsertMappingPath(barcode, path);
                    }
                    counter++;
                    VERBOSE_POWER_T2(counter, 100, "Processed " << counter << " reads.");
                }
            }
            INFO("FillMap finished")
            //INFO("Number of barcodes: " + std::to_string(barcode_codes_.GetSize()))
        }

        void FillMap(const lib_vector_t& libs_10x) {
            InitialFillMap();
            FillMapFrom10XReads(libs_10x);
            return;
        }

        virtual void InitialFillMap() = 0;

    protected:

        void InsertEntry(const EdgeId &edge, barcode_entry_t &entry) {
            auto key_and_value = std::make_pair(edge, entry);
            mapper_->edge_to_entry_.insert({edge, entry});
        }


        //todo works with a certain format of 10x only
        string GetTenXBarcodeFromRead(const io::SingleRead &read) {
            const size_t barcode_len = 16;
            const string prefix_string = "BC:Z";
            size_t prefix_len = prefix_string.size();
            size_t start_pos = read.name().find(prefix_string);
            if (start_pos != string::npos) {
                VERIFY(start_pos + prefix_len + barcode_len <= read.name().length())
                string barcode = read.name().substr(start_pos + 5, barcode_len);
                return barcode;
            }
            return "";
        }

        KmerMultiset BuildKmerMultisetFromStream(shared_ptr <io::ReadStream<io::PairedRead>> read_stream) {
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

        void ExtractKmersFromSeq(KmerMultiset &kmer_multiset,
                                 const io::SingleRead &read, size_t kmer_len) {
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
            mapper_->edge_to_entry_.at(edge).InsertBarcode(barcode, count, range);
        }

        bool IsAtEdgeTail(const EdgeId &edge, const omnigraph::Range &range) {
            return range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool IsAtEdgeHead(const omnigraph::Range &range) {
            return range.end_pos < tail_threshold_;
        }

        void InsertMappingPath(const BarcodeId &barcode, const MappingPath <EdgeId> &path) {
            for (size_t j = 0; j < path.size(); j++) {
                InsertBarcodeWithRange(barcode, path[j].first, path[j].second.mapped_range, 1);
            }
        }

        void InsertBarcodeWithRange(const BarcodeId &barcode, const EdgeId &edge,
                                    const omnigraph::Range &range, size_t count) {
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

        vector <io::SingleStreamPtr> GetStreamsFromLibs(const lib_vector_t& libs_10x) { ;
            vector <io::SingleStreamPtr> result;
            for (const auto& lib: libs_10x) {
                VERIFY(lib.type() == io::LibraryType::Clouds10x);
                for (auto it = lib.reads_begin(); it != lib.reads_end(); ++it) {
                    auto stream = io::EasyStream(*it, false);
                    result.push_back(stream);
                }
            }
            return result;
        }

        vector <vector<tslr_barcode_library>> SplitLibrary(const vector <tslr_barcode_library> &lib_vec,
                                                           size_t bucket_size) {
            size_t lib_size = lib_vec.size();
            size_t buckets = lib_size / bucket_size;
            vector <vector<tslr_barcode_library>> bucket_vec(buckets);
            for (size_t i = 0; i < buckets; ++i) {
                for (size_t j = 0; j < bucket_size; ++j) {
                    bucket_vec[i].push_back(lib_vec[i * bucket_size + j]);
                }
            }
            vector <tslr_barcode_library> last_bucket;
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
    };

    class SimpleMapperBuilder : public BarcodeIndexBuilder<SimpleEdgeEntry> {
        using BarcodeIndexBuilder::g_;
        using BarcodeIndexBuilder::mapper_;
    public:
        SimpleMapperBuilder(const Graph &g, const size_t tail_threshold) :
                BarcodeIndexBuilder(g, tail_threshold) {}

        void InitialFillMap() override {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                SimpleEdgeEntry entry(*it);

                InsertEntry(*it, entry);
            }
        }
    };

    class FrameMapperBuilder : public BarcodeIndexBuilder<FrameEdgeEntry> {
        using BarcodeIndexBuilder::g_;
        using BarcodeIndexBuilder::mapper_;
        const size_t frame_size_;
    public:
        FrameMapperBuilder(const Graph &g, const size_t tail_threshold, const size_t frame_size) :
                BarcodeIndexBuilder(g, tail_threshold),
                frame_size_(frame_size) {}

        void InitialFillMap() override {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                FrameEdgeEntry entry(*it, g_.length(*it), frame_size_);
                InsertEntry(*it, entry);
            }
        }
    };
}
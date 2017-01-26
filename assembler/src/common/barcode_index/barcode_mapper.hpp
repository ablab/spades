#pragma once

#include <boost/unordered_map.hpp>
#include <memory>
#include <utility>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <bitset>
#include "io/reads/paired_readers.hpp"
#include <common/assembly_graph/paths/mapping_path.hpp>
#include <cassert>
#include "common/modules/alignment/edge_index.hpp"
#include "common/modules/alignment/kmer_mapper.hpp"
#include "common/modules/alignment/sequence_mapper.hpp"
#include "common/pipeline/config_struct.hpp"
#include "common/utils/indices/edge_index_builders.hpp"
#include "common/utils/range.hpp"

using std::string;
using std::istringstream;
using namespace omnigraph;

namespace tslr_resolver {
    //constexpr int16_t max_barcodes = 384;

    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef debruijn_graph::EdgeIndex<Graph> Index;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    typedef debruijn_graph::KmerMapper<Graph> KmerSubs;
    typedef string BarcodeId;
    typedef RtSeq Kmer;
    typedef typename debruijn_graph::KmerFreeEdgeIndex<Graph, debruijn_graph::DefaultStoring> InnerIndex;
    typedef typename InnerIndex::KeyWithHash KeyWithHash;
    typedef typename debruijn_graph::EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;

    enum BarcodeLibraryType {
        TSLR,
        TenX,
        Unknown
    };

    inline BarcodeLibraryType GetLibType(const string type) {
        if (type == "tslr")
            return BarcodeLibraryType::TSLR;
        if (type == "tenx")
            return BarcodeLibraryType::TenX;
        return BarcodeLibraryType::Unknown;
    }

    struct tslr_barcode_library {
        string left_;
        string right_;
        string barcode_;
    };


    class BarcodeEncoder {
        std::unordered_map <BarcodeId, int64_t> codes_;
        int64_t barcode_encoder_size;
    public:
        BarcodeEncoder() :
                codes_(), barcode_encoder_size (0)
        { }

        void AddBarcode(const string &barcode) {
            auto it = codes_.find(barcode);
            if (it == codes_.end()) {
                codes_[barcode] = barcode_encoder_size;
                barcode_encoder_size++;
            }
        }

        int64_t GetCode (const string& barcode) const {
            VERIFY(codes_.find(barcode) != codes_.end());
            return codes_.at(barcode);
        }

        int64_t GetSize() const {
            return barcode_encoder_size;
        }
    };

    class KmerMultiset {
        boost::unordered_map<Kmer, size_t, Kmer::hash> storage_;

    public:
        KmerMultiset() : storage_() {}

        void Insert(const Kmer& kmer) {
            if (kmer.IsMinimal()) {
                if (storage_.find(kmer) == storage_.end()) {
                    storage_[kmer] = 1;
                }
                else {
                    ++storage_[kmer];
                }
            }
        }

        size_t size() const {
            return storage_.size();
        }

        decltype(storage_.cbegin()) cbegin() const {
            return storage_.cbegin();
        }

        decltype(storage_.cend()) cend() const {
            return storage_.cend();
        }
    };

    template <class barcode_entry_t>
    class HeadTailMapperBuilder;

    /*This structure contains barcode multiset extracted from reads aligned to the
     * beginning of edges in the assembly graph. */
    class BarcodeMapper {
    public:
    protected:
        const Graph& g_;
        size_t barcodes_number_ = 0;
    public:
        BarcodeMapper (const Graph &g) :
                g_(g), barcodes_number_(0) {}
        virtual ~BarcodeMapper() {}

        virtual size_t GetNumberOfBarcodes() const = 0;

        //Number of entries in the barcode map. Currently equals to number of edges.
        virtual size_t size() const = 0;

        virtual vector<int64_t> GetIntersection(const EdgeId& edge1, const EdgeId& edge2) const = 0;

        //Get number of shared barcodes between two edges.
        virtual size_t GetIntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double GetIntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double GetIntersectionSizeNormalizedByFirst(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual size_t GetUnionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;

        //Average barcode coverage of long edges
        virtual double AverageBarcodeCoverage () const = 0;

        //Number of barcodes on the beginning/end of the edge
        virtual size_t GetHeadBarcodeNumber(const EdgeId& edge) const = 0;
        virtual size_t GetTailBarcodeNumber(const EdgeId& edge) const = 0;

        //fixme these methods should be moved to DataScanner
        virtual void ReadEntry(ifstream& fin, const EdgeId& edge) = 0;
        virtual void WriteEntry(ofstream& fin, const EdgeId& edge) = 0;

        //Remove low abundant barcodes
        virtual void Filter(size_t abundancy_threshold, size_t gap_threshold) = 0;

        //Serialize barcode abundancies. Format:
        //abundancy: number of edges.
        virtual void SerializeOverallDistribution(const string& path) const = 0;
        virtual bool IsEmpty() = 0;

    };

    template <class barcode_entry_t>
    class HeadTailBarcodeMapper : public BarcodeMapper {
    friend class HeadTailMapperBuilder<barcode_entry_t>;
    //temporary?
    friend class BarcodeStatisticsCollector;
    protected:
        typedef std::unordered_map <EdgeId, barcode_entry_t> barcode_map_t;
        using BarcodeMapper::g_;
        using BarcodeMapper::barcodes_number_;
        barcode_map_t edge_to_distribution_;

    public:
        HeadTailBarcodeMapper (const Graph &g) :
                BarcodeMapper(g),
                edge_to_distribution_()
        {
            InitialFillMap();
        }

        HeadTailBarcodeMapper (const HeadTailBarcodeMapper& other) = default;

        virtual ~HeadTailBarcodeMapper() {}

        void InitialFillMap() {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_entry_t set(*it);
                edge_to_distribution_.insert({*it, set});
            }
        }

        size_t GetNumberOfBarcodes() const override {
            return barcodes_number_;
        }

        size_t size() const {
            return edge_to_distribution_.size();
        }

        typename barcode_map_t::const_iterator cbegin() const noexcept {
            return edge_to_distribution_.cbegin();
        }

        typename barcode_map_t::const_iterator cend() const noexcept {
            return edge_to_distribution_.cend();
        }


        virtual double GetIntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const override {
            if (GetUnionSize(edge1, edge2)) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetUnionSize(edge1, edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId& edge1, const EdgeId& edge2) const override {
            if (GetHeadBarcodeNumber(edge2) > 0) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetHeadBarcodeNumber(edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedByFirst(const EdgeId& edge1, const EdgeId& edge2) const override {
            if (GetTailBarcodeNumber(edge1) > 0) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetTailBarcodeNumber(edge1));
            }
            return 0;
        }


        size_t GetHeadBarcodeNumber(const EdgeId& edge) const override {
            return GetEntryHeads(edge).Size();
        }

        size_t GetTailBarcodeNumber(const EdgeId& edge) const override {
            return GetEntryTails(edge).Size();
        }

        bool IsEmpty() override {
            const double empty_threshold = 0.00001;
            return AverageBarcodeCoverage() < empty_threshold;
        }

        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall = 0;
            int64_t long_edges = 0;
            size_t len_threshold = cfg::get().ts_res.len_threshold;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > len_threshold) {
                    long_edges++;
                    barcodes_overall += GetTailBarcodeNumber(*it);
                }
            }
            DEBUG("tails: " + std::to_string(barcodes_overall));
            DEBUG("Long edges" + long_edges);
            return static_cast <double> (barcodes_overall) / static_cast <double> (long_edges);
        }

        vector<int64_t> GetIntersection(const EdgeId& edge1, const EdgeId& edge2) const {
            return GetEntryTails(edge1).GetIntersection(GetEntryHeads(edge2));
        }

        size_t GetIntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const override {
            return GetEntryTails(edge1).GetIntersectionSize(GetEntryHeads(edge2));
        }

        size_t GetUnionSize(const EdgeId& edge1, const EdgeId& edge2) const override {
            return GetEntryTails(edge1).GetUnionSize(GetEntryHeads(edge2));
        }


        //Delete low abundant barcodes from every edge
        void Filter(size_t trimming_threshold, size_t gap_threshold) override {
            for (auto entry = edge_to_distribution_.begin(); entry != edge_to_distribution_.end(); ++entry) {
                entry->second.Filter(trimming_threshold, gap_threshold);
            }
        }

        void SerializeOverallDistribution(const string& path) const override {
            ofstream fout;
            fout.open(path);
            std::map <size_t, size_t> overall_distr;
            INFO("Serializing distribution")
            for (const auto& entry: edge_to_distribution_) {
                //fixme config
                if (g_.length(entry.first) > cfg::get().ts_res.len_threshold) {
                    const auto &current_distr = edge_to_distribution_.at(entry.first);
                    for (auto it = current_distr.cbegin();
                         it != current_distr.cend(); ++it) {
                        overall_distr[it->second.GetCount()]++;
                    }
                }
            }
            for (const auto& entry : overall_distr) {
                fout << entry.first << ": " << entry.second << endl;
            }
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge) override {
            barcode_entry_t entry(edge);
            entry.Deserialize(fin);
            edge_to_distribution_[edge] = entry;
            DEBUG(edge.int_id());
            DEBUG(entry.Size());
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge) override {
            fout << g_.int_id(edge) << std::endl;
            GetEntryHeads(edge).Serialize(fout);
        }

        barcode_entry_t GetEntryHeads(const EdgeId& edge) const {
            return edge_to_distribution_.at(edge);
        }

        barcode_entry_t GetEntryTails(const EdgeId& edge) const {
            return edge_to_distribution_.at(g_.conjugate(edge));
        }

    };

    class BarcodeInfo {
        size_t count_;
        Range range_;
    public:
        BarcodeInfo(): count_(0), range_() {}
        BarcodeInfo(size_t count, const Range& range): count_(count), range_(range) {}

        bool IsFarFromEdgeHead(size_t gap_threshold) {
            return range_.start_pos > gap_threshold;
        }

        bool IsLowReadCount(size_t trimming_threshold) {
            return count_ < trimming_threshold;
        }

        void Update(const BarcodeInfo& info) {
            count_ += info.count_;
            range_.start_pos = std::min(range_.start_pos, info.range_.start_pos);
            range_.end_pos = std::max(range_.end_pos, info.range_.end_pos);
        }

        size_t GetCount() const {
            return count_;
        }


        friend ostream& operator <<(ostream& os, const BarcodeInfo& info);
        friend istream& operator >>(istream& is, BarcodeInfo& info);
    };

    inline ostream& operator <<(ostream& os, const BarcodeInfo& info)
    {
        os << info.count_ << " " << info.range_.start_pos << " " << info.range_.end_pos;
        return os;
    }

    inline istream& operator >>(istream& os, BarcodeInfo& info)
    {
        size_t count;
        size_t range_start;
        size_t range_end;
        Range range;
        os >> count;
        os >> range_start;
        os >> range_end;
        info.count_ = count;
        info.range_ = Range(range_start, range_end);
        return os;
    }

    //Contains abundancy for each barcode aligned to given edge
    class EdgeEntry {
    protected:
        friend class HeadTailBarcodeMapper<EdgeEntry>;
        friend class HeadTailMapperBuilder<EdgeEntry>;

        typedef std::unordered_map <int64_t, BarcodeInfo> barcode_distribution_t;
        EdgeId edge_;
        barcode_distribution_t barcode_distribution_;

    public:
        EdgeEntry():
            edge_(), barcode_distribution_() {};
        EdgeEntry(const EdgeId& edge) :
                edge_(edge), barcode_distribution_() {}

        virtual ~EdgeEntry() {}


        barcode_distribution_t GetDistribution() const {
            return barcode_distribution_;
        }

        EdgeId GetEdge() const {
            return edge_;
        }

        vector <int64_t> GetIntersection(const EdgeEntry& other) const {
            vector <int64_t> result;
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end(); ++it) {
                if (other.GetDistribution().find(it-> first) != other.GetDistribution().end()) {
                    result.push_back(it->first);
                }
            }
            return result;
        }

        size_t GetIntersectionSize(const EdgeEntry &other) const {
            return GetIntersection(other).size();
        }

        size_t GetUnionSize(const EdgeEntry& other) const {
            auto distr_this = barcode_distribution_;
            auto distr_other = other.GetDistribution();
            return Size() + other.Size() - GetIntersectionSize(other);
        }

        void InsertSet (const barcode_distribution_t& set) {
            barcode_distribution_ = set;
        }

        void Filter(size_t trimming_threshold, size_t gap_threshold) {
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end() ;) {
                if (it->second.IsLowReadCount(trimming_threshold) or
                        it->second.IsFarFromEdgeHead(gap_threshold)) {
                    barcode_distribution_.erase(it++);
                }
                else {
                    ++it;
                }
            }
        }

        size_t Size() const {
            return barcode_distribution_.size();
        }

        virtual void Serialize(ofstream& fout) {
            SerializeDistribution(fout, barcode_distribution_);
        }

        virtual void Deserialize(ifstream& fin) {
            DeserializeDistribution(fin, barcode_distribution_);
        }

        decltype(barcode_distribution_.cbegin()) cbegin() const {
            return barcode_distribution_.cbegin();
        }

        decltype(barcode_distribution_.cend()) cend() const {
            return barcode_distribution_.cend();
        }

    protected:
        void SerializeDistribution(ofstream &fout, const barcode_distribution_t &distribution) {
            //INFO("Serializing entry")
            fout << distribution.size() << endl;
            for (auto entry : distribution) {
                fout << entry.first << ' ' << entry.second << endl;
            }
        }

        void DeserializeDistribution(ifstream &fin, barcode_distribution_t &distribution) {
            //INFO("Deserializing entry")
            size_t distr_size;
            fin >> distr_size;
            //INFO(distr_size)
            for (size_t i = 0; i < distr_size; ++i) {
                int64_t bid;
                BarcodeInfo info;
                fin >> bid >> info;
                InsertBarcode(distribution, bid, info);
            }
        }

    private:


        void InsertBarcode(int64_t code, const BarcodeInfo& info) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                barcode_distribution_.insert({code, info});
            }
            else {
                barcode_distribution_.at(code).Update(info);
            }
        }

        void InsertBarcode(barcode_distribution_t& distribution, int64_t code, const BarcodeInfo& info) {
            if (distribution.find(code) == distribution.end()) {
                distribution.insert({code, info});
            }
            else {
                distribution.at(code).Update(info);
            }
        }
    };


    template <class barcode_entry_t>
    class HeadTailMapperBuilder {
        const Graph& g_;
        shared_ptr<HeadTailBarcodeMapper<barcode_entry_t>> mapper_;
        size_t tail_threshold_;
        BarcodeEncoder barcode_codes_;

    public:
        HeadTailMapperBuilder(const Graph& g, size_t tail_threshold) :
                g_(g),
                mapper_(make_shared<HeadTailBarcodeMapper<barcode_entry_t>>(g)),
                tail_threshold_(tail_threshold),
                barcode_codes_() {}
        ~HeadTailMapperBuilder() {}

        DECL_LOGGER("BarcodeMapperBuilder")

        shared_ptr<HeadTailBarcodeMapper<barcode_entry_t>> GetMapper() {
            return mapper_;
        }

        void FillMapFromDemultiplexedDataset(const Index &index, const KmerSubs &kmer_mapper) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            auto lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);

            //Process every barcode from truspades dataset
            for (size_t i = 0; i < lib_vec.size(); ++i) {
                std::string barcode = lib_vec[i].barcode_;
                std::shared_ptr<io::ReadStream<io::PairedRead>> paired_stream =
                        make_shared<io::SeparatePairedReadStream> (lib_vec[i].left_, lib_vec[i].right_, 1);
                io::PairedRead read;
                while (!paired_stream->eof()) {
                    *paired_stream >> read;
                    auto path_first = mapper -> MapRead(read.first());
                    auto path_second = mapper -> MapRead(read.second());
                    InsertMappingPath(barcode, path_first);
                    InsertMappingPath(barcode, path_second);
                }
                VERBOSE_POWER_T2(i, 100,
                                 "Processed " << i << " barcodes from " << lib_vec.size() << " (" << i * 100 / lib_vec.size()
                                              << "%)");
                if (lib_vec.size() > 10 && i % (lib_vec.size() / 10 + 1) == 0) {
                    INFO("Processed " << i << " barcodes from " << lib_vec.size() << " (" << i * 100 / lib_vec.size() << "%)");
                }
                ++mapper_->barcodes_number_;
            }
        }

        void FillMapUsingKmerMultisetParallel(const Index& index, const KmerSubs& kmer_mapper, size_t n_threads) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            const auto& lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);
            const auto& bucket_vec = SplitLibrary(lib_vec, n_threads);
            INFO("Library splitted");
            for (size_t i = 0; i < bucket_vec.size(); ++i) {
#pragma omp parallel for num_threads(n_threads)
                for (size_t j = 0; j < bucket_vec[i].size(); ++j) {
                    const auto& lib = bucket_vec[i][j];
                    std::string barcode = lib.barcode_;

                    std::shared_ptr<io::ReadStream<io::PairedRead>> paired_stream =
                            make_shared<io::SeparatePairedReadStream>(lib.left_, lib.right_, 1);
                    const KmerMultiset& kmer_multiset = BuildKmerMultisetFromStream(paired_stream);

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

        void FillMapUsingSubIndex (const Index& index, const KmerSubs& kmer_mapper) {
            //fixme move to command line
            std::string tslr_dataset = cfg::get().ts_res.read_cloud_dataset;

            auto lib_vec = GetLibrary(tslr_dataset);
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);
            size_t global_counter = 0;
            size_t counter = 0;

            //Process every barcode from truspades dataset
            for (size_t i = 0; i < lib_vec.size(); ++i) {
                std::string barcode = lib_vec[i].barcode_;
                Index barcode_subindex(g_, cfg::get().tmp_dir);
                InnerIndex subindex = InnerIndex(g_, cfg::get().tmp_dir);

                io::ReadStreamList<io::SingleRead> streams;
                streams.push_back(io::EasyStream(lib_vec[i].left_, false /*followed_by_rc*/));
                streams.push_back(io::EasyStream(lib_vec[i].right_, false /*followed_by_rc*/));
                IndexBuilder().BuildIndexFromStream(subindex, streams, 0);

                INFO("Extracting kmer multiset using subindex")

                for (auto it = subindex.kmer_begin(); it.good(); ++it) {
                    Kmer kmer = Kmer(g_.k() + 1, *it);
                    KeyWithHash kh = subindex.ConstructKWH(kmer);
                    const auto& edge_and_offset = index.get(kmer);
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

        void FillMapFrom10XReads(const Index &index, const KmerSubs &kmer_mapper) {
            INFO("Starting barcode index construction from 10X reads")
            std::string read_cloud_dataset = cfg::get().ts_res.read_cloud_dataset;
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);

            auto streams = GetStreamsFromFile(read_cloud_dataset);
            //Process every read from 10X dataset
            io::SingleRead read;
            size_t counter = 0;
            for (auto stream: streams) {
                while (!stream->eof()) {
                    *stream >> read;
                    auto barcode = GetBarcodeFromRead(read);
                    if (barcode != "") {
                        barcode_codes_.AddBarcode(barcode);
                        const auto& path = mapper->MapRead(read);
                        //INFO("Inserting")
//                        INFO(path.size())
                        InsertMappingPath(barcode, path);
                    }
                    counter++;
                    VERBOSE_POWER_T2(counter, 100, "Processed " << counter << " reads.");
                }
            }
            INFO("FillMap finished")
            //INFO("Number of barcodes: " + std::to_string(barcode_codes_.GetSize()))
        }

        void FillMap(BarcodeLibraryType lib_type, const Index& index, const KmerSubs& kmer_mapper, size_t nthreads) {
            switch(lib_type) {
                case TSLR :
                    //FillMapUsingKmerMultisetParallel(index, kmer_mapper, nthreads);
                    FillMapFromDemultiplexedDataset(index, kmer_mapper);
                    break;
                case TenX :
                    FillMapFrom10XReads(index, kmer_mapper);
                    break;
                default:
                    WARN("Unknown library type, failed to fill barcode map.");
                    return;
            }
        }

    protected:

        //todo works with a certain format of 10x only
        std::string GetBarcodeFromRead(const io::SingleRead &read) {
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

        KmerMultiset BuildKmerMultisetFromStream(shared_ptr<io::ReadStream<io::PairedRead>> read_stream) {
            KmerMultiset kmer_multiset;
            size_t read_counter = 0;
            io::PairedRead read;
            size_t kmer_len = g_.k() + 1;
            while (!read_stream -> eof()) {
                read_counter += 2;
                *read_stream >> read;
                ExtractKmersFromSeq(kmer_multiset, read.first(), kmer_len);
                ExtractKmersFromSeq(kmer_multiset, read.second(), kmer_len);
            }
            return kmer_multiset;
        }

        void ExtractKmersFromSeq(KmerMultiset& kmer_multiset,
                          const io::SingleRead& read, size_t kmer_len) {
            if (read.IsValid()) {
                const Sequence& seq = read.sequence();
                Kmer kmer = seq.start<Kmer>(kmer_len);
                for (size_t j = kmer_len; j < seq.size(); ++j) {
                    kmer_multiset.Insert(kmer);
                    kmer <<= seq[j];
                }
                kmer_multiset.Insert(kmer);
            }
        }


        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge, const BarcodeInfo& info) {
            int64_t code = barcode_codes_.GetCode(barcode);
            mapper_ -> edge_to_distribution_.at(edge).InsertBarcode(code, info);
        }

        bool IsAtEdgeTail(const EdgeId& edge, const omnigraph::Range &range) {
            return range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool IsAtEdgeHead(const omnigraph::Range &range) {
            return range.end_pos < tail_threshold_;
        }

        void InsertMappingPath(const BarcodeId& barcode, const MappingPath<EdgeId>& path) {
            for (size_t j = 0; j < path.size(); j++) {
                InsertBarcodeWithRange(barcode, path[j].first, path[j].second.mapped_range, 1);
            }
        }

        void InsertBarcodeWithRange(const BarcodeId& barcode, const EdgeId& edge,
                                    const omnigraph::Range& range, size_t count) {
            BarcodeInfo info(count, range);
            if (IsAtEdgeHead(range))
                InsertBarcode(barcode, edge, info);
            if (IsAtEdgeTail(edge, range))
                InsertBarcode(barcode, g_.conjugate(edge), info);
        }


        std::vector <tslr_barcode_library> GetLibrary(const string& tslr_dataset_name) {
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

        vector <io::SingleStreamPtr> GetStreamsFromFile(const string& filename) {
            std::ifstream fin;
            fin.open(filename);
            string read_filename;
            vector <io::SingleStreamPtr> result;
            while(getline(fin, read_filename)) {
                auto stream = io::EasyStream(read_filename, false);
                result.push_back(stream);
            }
            return result;
        }

        vector <vector <tslr_barcode_library>> SplitLibrary(const vector <tslr_barcode_library>& lib_vec,
                                                            size_t bucket_size) {
            size_t lib_size = lib_vec.size();
            size_t buckets = lib_size / bucket_size;
            vector<vector <tslr_barcode_library>> bucket_vec(buckets);
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
            for (const auto& vec : bucket_vec) {
                sum_size += vec.size();
            }
            VERIFY(lib_size == sum_size);
            return bucket_vec;
        }
    };

} //tslr_resolver

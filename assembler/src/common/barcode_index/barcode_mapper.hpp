#pragma once

#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>
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
#include <common/modules/alignment/bwa_sequence_mapper.hpp>
#include "common/modules/alignment/edge_index.hpp"
#include "common/modules/alignment/kmer_mapper.hpp"
#include "common/modules/alignment/sequence_mapper.hpp"
#include "common/pipeline/config_struct.hpp"
#include "common/utils/indices/edge_index_builders.hpp"
#include "common/utils/range.hpp"

using std::string;
using std::istringstream;
using namespace omnigraph;

namespace barcode_index {
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
    class BarcodeIndexBuilder;

    template <class barcode_entry_t>
    class BarcodeIndexInfoExtractor;

    /*This structure contains barcode multiset extracted from reads aligned to the
     * beginning of edges in the assembly graph. */
    class AbstractBarcodeIndex {
    public:
    protected:
        const Graph& g_;
        size_t barcodes_number_ = 0;
    public:
        AbstractBarcodeIndex (const Graph &g) :
                g_(g), barcodes_number_(0) {}
        virtual ~AbstractBarcodeIndex() {}

        virtual size_t GetNumberOfBarcodes() const = 0;

        //Number of entries in the barcode map. Currently equals to number of edges.
        virtual size_t size() const = 0;

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
        virtual bool IsEmpty() = 0;

    };

    template <class barcode_entry_t>
    class BarcodeIndex : public AbstractBarcodeIndex {
    friend class BarcodeIndexBuilder<barcode_entry_t>;
    friend class BarcodeIndexInfoExtractor<barcode_entry_t>;
    friend class BarcodeStatisticsCollector;
    protected:
        typedef std::unordered_map <EdgeId, barcode_entry_t> barcode_map_t;
        using AbstractBarcodeIndex::g_;
        using AbstractBarcodeIndex::barcodes_number_;
        barcode_map_t edge_to_entry_;

    public:
        BarcodeIndex (const Graph &g) :
                AbstractBarcodeIndex(g),
                edge_to_entry_()
        {}

        BarcodeIndex (const BarcodeIndex& other) = default;

        virtual ~BarcodeIndex() {}

        void InitialFillMap() {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_entry_t set(*it);
                edge_to_entry_.insert({*it, set});
            }
        }

        size_t GetNumberOfBarcodes() const override {
            return barcodes_number_;
        }

        size_t size() const {
            return edge_to_entry_.size();
        }

        typename barcode_map_t::const_iterator cbegin() const noexcept {
            return edge_to_entry_.cbegin();
        }

        typename barcode_map_t::const_iterator cend() const noexcept {
            return edge_to_entry_.cend();
        }


        size_t GetHeadBarcodeNumber(const EdgeId& edge) const override {
            return GetEntryHeads(edge).Size();
        }

        size_t GetTailBarcodeNumber(const EdgeId& edge) const override {
            return GetEntryTails(edge).Size();
        }

        bool IsEmpty() override {
            return size() == 0;
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
            INFO("Barcodes: " << barcodes_overall);
            return static_cast <double> (barcodes_overall) / static_cast <double> (long_edges);
        }

        //Delete low abundant barcodes from every edge
        void Filter(size_t trimming_threshold, size_t gap_threshold) override {
            for (auto entry = edge_to_entry_.begin(); entry != edge_to_entry_.end(); ++entry) {
                entry->second.Filter(trimming_threshold, gap_threshold);
            }
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge) override {
            edge_to_entry_[edge].Deserialize(fin);
            DEBUG(edge.int_id());
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge) override {
            fout << g_.int_id(edge) << std::endl;
            GetEntryHeads(edge).Serialize(fout);
        }

        typename barcode_map_t::const_iterator GetEntryTailsIterator(const EdgeId& edge) const {
            return edge_to_entry_.find(g_.conjugate(edge));
        }

        typename barcode_map_t::const_iterator GetEntryHeadsIterator(const EdgeId& edge) const {
            return edge_to_entry_.find(edge);
        }

        barcode_entry_t GetEntryHeads(const EdgeId& edge) const {
            return edge_to_entry_.at(edge);
        }

        barcode_entry_t GetEntryTails(const EdgeId& edge) const {
            return edge_to_entry_.at(g_.conjugate(edge));
        }
    };

    class SimpleBarcodeInfo {
        size_t count_;
        Range range_;
    public:
        SimpleBarcodeInfo(): count_(0), range_() {}
        SimpleBarcodeInfo(size_t count, const Range& range): count_(count), range_(range) {}

        void Update(size_t count, const Range& range) {
            count_ += count;
            range_.start_pos = std::min(range_.start_pos, range.start_pos);
            range_.end_pos = std::max(range_.end_pos, range.end_pos);
        }

        void Update(const SimpleBarcodeInfo& other) {
            count_ += other.GetCount();
            Range range;
            range_.start_pos = std::min(range_.start_pos, other.GetRange().start_pos);
            range_.end_pos = std::max(range_.end_pos, other.GetRange().end_pos);
        }

        size_t GetCount() const {
            return count_;
        }

        Range GetRange() const {
            return range_;
        }
        friend ostream& operator <<(ostream& os, const SimpleBarcodeInfo& info);
        friend istream& operator >>(istream& is, SimpleBarcodeInfo& info);
    };

    inline ostream& operator <<(ostream& os, const SimpleBarcodeInfo& info)
    {
        os << info.count_ << " " << info.range_.start_pos << " " << info.range_.end_pos;
        return os;
    }

    inline istream& operator >>(istream& os, SimpleBarcodeInfo& info)
    {
        size_t range_start;
        size_t range_end;
        os >> info.count_;
        os >> range_start;
        os >> range_end;
        info.range_ = Range(range_start, range_end);
        return os;
    }

    class FrameBarcodeInfo {
        size_t count_;
        boost::dynamic_bitset<> is_on_frame_;
        size_t leftmost_index_;
        size_t rightmost_index_;
    public:

        FrameBarcodeInfo(size_t frames = 0): count_(0), is_on_frame_(), leftmost_index_(frames), rightmost_index_(0) {
            is_on_frame_.resize(frames, false);
        }

        void Update(size_t count, size_t left_frame, size_t right_frame) {
            count_ += count;
            for (size_t i = left_frame; i <= right_frame; ++i) {
                is_on_frame_.set(i);
            }
            leftmost_index_ = std::min(left_frame, leftmost_index_);
            rightmost_index_ = std::max(right_frame, rightmost_index_);
        }

        void Update(const FrameBarcodeInfo& other) {
            is_on_frame_ |= other.is_on_frame_;
            leftmost_index_ = std::min(leftmost_index_, other.leftmost_index_);
            rightmost_index_ = std::max(leftmost_index_, other.rightmost_index_);
            count_ += other.count_;
        }

        size_t GetCount() const {
            return count_;
        }

        size_t GetLeftMost() const {
            return leftmost_index_;
        }

        size_t GetRightMost() const {
            return rightmost_index_;
        }

        bool GetFrame(size_t frame) const {
            return is_on_frame_[frame];
        }

        friend ostream& operator <<(ostream& os, const FrameBarcodeInfo& info);
        friend istream& operator >>(istream& is, FrameBarcodeInfo& info);
    };

    inline ostream& operator <<(ostream& os, const FrameBarcodeInfo& info)
    {
        os << info.count_ << " " << info.is_on_frame_;
        return os;
    }

    inline istream& operator >>(istream& os, FrameBarcodeInfo& info)
    {
        os >> info.count_;
        os >> info.is_on_frame_;
        info.leftmost_index_ = info.is_on_frame_.find_first();
        size_t rightmost = 0;
        for (size_t i = info.is_on_frame_.size() - 1; i > 0; --i) {
            if (info.is_on_frame_.test(i)) {
                rightmost = i;
                break;
            }
        }
        info.rightmost_index_ = rightmost;
        return os;
    }

    template <class barcode_info_t>
    class EdgeEntry {
    protected:
        typedef std::unordered_map <int64_t, barcode_info_t> barcode_distribution_t;
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

        //fixme move to info extractor
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
            size_t result = 0;
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end(); ++it) {
                if (other.GetDistribution().find(it-> first) != other.GetDistribution().end()) {
                    result++;
                }
            }
            return result;
        }

        size_t GetUnionSize(const EdgeEntry& other) const {
            auto distr_this = barcode_distribution_;
            auto distr_other = other.GetDistribution();
            return Size() + other.Size() - GetIntersectionSize(other);
        }

        void InsertSet (const barcode_distribution_t& set) {
            barcode_distribution_ = set;
        }

        size_t Size() const {
            return barcode_distribution_.size();
        }

        virtual void Serialize(ofstream& fout) {
            SerializeDistribution(fout);
        }

        virtual void Deserialize(ifstream& fin) {
            DeserializeDistribution(fin);
        }

        typename barcode_distribution_t::const_iterator begin() const {
            return barcode_distribution_.cbegin();
        }

        typename barcode_distribution_t::const_iterator end() const {
            return barcode_distribution_.cend();
        }

        typename barcode_distribution_t::const_iterator cbegin() const {
            return barcode_distribution_.cbegin();
        }

        typename barcode_distribution_t::const_iterator cend() const {
            return barcode_distribution_.cend();
        }

        bool has_barcode(int64_t barcode) const {
            return barcode_distribution_.find(barcode) != barcode_distribution_.end();
        }

        typename barcode_distribution_t::const_iterator get_barcode(int64_t barcode) const {
            return barcode_distribution_.find(barcode);
        }

    protected:
        void SerializeDistribution(ofstream &fout) {
            //INFO("Serializing entry")
            fout << barcode_distribution_.size() << endl;
            for (auto entry : barcode_distribution_) {
                fout << entry.first << ' ' << entry.second << endl;
            }
        }

        void DeserializeDistribution(ifstream &fin) {
            //INFO("Deserializing entry")
            size_t distr_size;
            fin >> distr_size;
            //INFO(distr_size)
            for (size_t i = 0; i < distr_size; ++i) {
                int64_t bid;
                barcode_info_t info;
                fin >> bid >> info;
                InsertInfo(bid, info);
            }
        }

        virtual void InsertInfo(int64_t code, const barcode_info_t& info) = 0;
        virtual void InsertBarcode(int64_t code, const size_t count, const Range& range) = 0;
    };

    //Contains abundancy for each barcode aligned to given edge
    class SimpleEdgeEntry : public EdgeEntry<SimpleBarcodeInfo> {
        friend class BarcodeIndex<SimpleEdgeEntry>;
        friend class BarcodeIndexBuilder<SimpleEdgeEntry>;
        friend class BarcodeIndexInfoExtractor<SimpleEdgeEntry>;
    protected:
        using EdgeEntry::barcode_distribution_t;
        using EdgeEntry::barcode_distribution_;
        using EdgeEntry::edge_;

    public:
        SimpleEdgeEntry():
            EdgeEntry() {}
        SimpleEdgeEntry(const EdgeId& edge) :
            EdgeEntry(edge) {}

        ~SimpleEdgeEntry() {}

        void Filter(size_t trimming_threshold, size_t gap_threshold) {
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end() ;) {
                if (IsLowReadCount(trimming_threshold, it->second) or
                        IsFarFromEdgeHead(gap_threshold, it->second)) {
                    barcode_distribution_.erase(it++);
                }
                else {
                    ++it;
                }
            }
        }

    protected:
        void InsertInfo(int64_t code, const SimpleBarcodeInfo &info) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                barcode_distribution_.insert({code, info});
            }
            else {
                barcode_distribution_.at(code).Update(info);
            }
        }

        void InsertBarcode(int64_t code, const size_t count, const Range& range) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                SimpleBarcodeInfo info(count, range);
                barcode_distribution_.insert({code, info});
            }
            else {
                barcode_distribution_.at(code).Update(count, range);
            }
        }


        bool IsFarFromEdgeHead(size_t gap_threshold, const SimpleBarcodeInfo& info) {
            return info.GetRange().start_pos > gap_threshold;
        }

        bool IsLowReadCount(size_t trimming_threshold, const SimpleBarcodeInfo& info) {
            return info.GetCount() < trimming_threshold;
        }
    };

    class FrameEdgeEntry : public EdgeEntry<FrameBarcodeInfo> {
        friend class BarcodeIndex<FrameEdgeEntry>;
        friend class BarcodeIndexBuilder<FrameEdgeEntry>;
        friend class BarcodeIndexInfoExtractor<FrameEdgeEntry>;
    protected:
        using EdgeEntry::barcode_distribution_t;
        using EdgeEntry::barcode_distribution_;
        using EdgeEntry::edge_;
        size_t edge_length_;
        size_t frame_size_;
        size_t number_of_frames_;

    public:
        FrameEdgeEntry():
            EdgeEntry(),
            edge_length_(0),
            frame_size_(0),
            number_of_frames_(0) {}
        FrameEdgeEntry(const EdgeId& edge, size_t edge_length, size_t frame_size) :
            EdgeEntry(edge),
            edge_length_(edge_length),
            frame_size_(frame_size),
            number_of_frames_(edge_length / frame_size + 1) {}

        ~FrameEdgeEntry() {}

        void Filter(size_t trimming_threshold, size_t gap_threshold) {
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end() ;) {
                if (IsLowReadCount(trimming_threshold, it->second) or
                        IsFarFromEdgeHead(gap_threshold, it->second)) {
                    barcode_distribution_.erase(it++);
                }
                else {
                    ++it;
                }
            }
        }

        size_t GetFrameSize() const {
            return frame_size_;
        }

    protected:
        void InsertInfo(int64_t code, const FrameBarcodeInfo &info) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                barcode_distribution_.insert({code, info});
            }
            else {
                barcode_distribution_.at(code).Update(info);
            }
        }

        void InsertBarcode(int64_t code, const size_t count, const Range& range) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                FrameBarcodeInfo info(number_of_frames_) ;
                barcode_distribution_.insert({code, info});
            }
            else {
                size_t left_frame = GetFrameFromPos(range.start_pos);
                size_t right_frame = GetFrameFromPos(range.end_pos);
                DEBUG("Range: " << range);
                DEBUG("Frames: " << left_frame << " " << right_frame);
                barcode_distribution_.at(code).Update(count, left_frame, right_frame);
            }
        }


        bool IsFarFromEdgeHead(size_t gap_threshold, const FrameBarcodeInfo& info) {
            return info.GetLeftMost() > gap_threshold / frame_size_;
        }

        bool IsLowReadCount(size_t trimming_threshold, const FrameBarcodeInfo& info) {
            return info.GetCount() < trimming_threshold;
        }

    private:
        //fixme last frame is larger than the others
        size_t GetFrameFromPos(size_t pos) {
            return pos / frame_size_;
        }

    };


    template <class barcode_entry_t>
    class BarcodeIndexBuilder {
    protected:
        const Graph& g_;
        shared_ptr<BarcodeIndex<barcode_entry_t>> mapper_;
        size_t tail_threshold_;
        BarcodeEncoder barcode_codes_;

    public:
        BarcodeIndexBuilder(const Graph& g, size_t tail_threshold) :
                g_(g),
                mapper_(make_shared<BarcodeIndex<barcode_entry_t>>(g)),
                tail_threshold_(tail_threshold),
                barcode_codes_() {}
        ~BarcodeIndexBuilder() {}

        DECL_LOGGER("BarcodeMapperBuilder")

        shared_ptr<BarcodeIndex<barcode_entry_t>> GetMapper() {
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

        void FillMapFrom10XReads() {
            INFO("Starting barcode index construction from 10X reads")
            std::string read_cloud_dataset = cfg::get().ts_res.read_cloud_dataset;
            INFO(read_cloud_dataset);
            auto mapper = std::make_shared<alignment::BWAReadMapper<Graph> >
                    (g_);

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
            InitialFillMap();
            switch(lib_type) {
                case TSLR :
                    //FillMapUsingKmerMultisetParallel(index, kmer_mapper, nthreads);
                    FillMapFromDemultiplexedDataset(index, kmer_mapper);
                    break;
                case TenX :
                    FillMapFrom10XReads();
                    break;
                default:
                    WARN("Unknown library type, failed to fill barcode map.");
                    return;
            }
        }

        virtual void InitialFillMap() = 0;

    protected:

        void InsertEntry (const EdgeId& edge, barcode_entry_t& entry) {
            auto key_and_value = std::make_pair(edge, entry);
            mapper_->edge_to_entry_.insert({edge, entry});
        }


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


        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge, size_t count, const Range& range) {
            int64_t code = barcode_codes_.GetCode(barcode);
            mapper_ -> edge_to_entry_.at(edge).InsertBarcode(code, count, range);
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
            if (IsAtEdgeHead(range))
                InsertBarcode(barcode, edge, count, range);
            if (IsAtEdgeTail(edge, range))
                InsertBarcode(barcode, g_.conjugate(edge), count, range.Invert(g_.length(edge)));
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

    class SimpleMapperBuilder : public BarcodeIndexBuilder<SimpleEdgeEntry> {
        using BarcodeIndexBuilder::g_;
        using BarcodeIndexBuilder::mapper_;
    public:
        SimpleMapperBuilder(const Graph& g, const size_t tail_threshold) :
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
        FrameMapperBuilder(const Graph& g, const size_t tail_threshold, const size_t frame_size) :
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

    class AbstractBarcodeIndexInfoExtractor {
    public:
        AbstractBarcodeIndexInfoExtractor() {}
        virtual ~AbstractBarcodeIndexInfoExtractor() {};
        virtual double GetIntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double GetIntersectionSizeNormalizedByFirst(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual size_t GetHeadBarcodeNumber(const EdgeId& edge) const = 0;
        virtual size_t GetTailBarcodeNumber(const EdgeId& edge) const = 0;
        virtual double AverageBarcodeCoverage() const = 0;
        virtual vector<int64_t> GetIntersection(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual size_t GetIntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual size_t GetUnionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual bool has_barcode(const EdgeId& edge, int64_t barcode) const = 0;
    };

    template <class barcode_entry_t>
    class BarcodeIndexInfoExtractor : public AbstractBarcodeIndexInfoExtractor {
    protected:
        shared_ptr<BarcodeIndex<barcode_entry_t>> mapper_;
        const Graph& g_;
    public:

        BarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper, const Graph& g):
        mapper_(std::dynamic_pointer_cast<BarcodeIndex<barcode_entry_t>>(abstract_mapper)),
        g_(g) {}

        double GetIntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const override {
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
            return mapper_->GetEntryHeads(edge).Size();
        }

        size_t GetTailBarcodeNumber(const EdgeId& edge) const override {
            return mapper_->GetEntryTails(edge).Size();
        }

        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall = 0;
            int64_t long_edges = 0;
            //fixme config
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

        vector<int64_t> GetIntersection(const EdgeId& edge1, const EdgeId& edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetIntersection(it_head->second);
        }

        size_t GetIntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetIntersectionSize(it_head->second);
        }

        size_t GetUnionSize(const EdgeId& edge1, const EdgeId& edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetUnionSize(it_head->second);
        }

        bool has_barcode(const EdgeId& edge, int64_t barcode) const override {
            return mapper_->GetEntryHeads(edge).has_barcode(barcode);
        }

        typename barcode_entry_t::barcode_distribution_t::const_iterator barcode_iterator_begin(const EdgeId& edge) {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            return entry_it->second.begin();
        }

        typename barcode_entry_t::barcode_distribution_t::const_iterator barcode_iterator_end(const EdgeId& edge) {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            return entry_it->second.end();
        }
    };

    class FrameBarcodeIndexInfoExtractor : public BarcodeIndexInfoExtractor<FrameEdgeEntry> {
    public:
        FrameBarcodeIndexInfoExtractor(shared_ptr <AbstractBarcodeIndex> abstract_mapper_ptr, const Graph& g) :
            BarcodeIndexInfoExtractor(abstract_mapper_ptr, g) {}

        size_t GetIntersectionSize(const EdgeId& first, const EdgeId& second, size_t gap_threshold) {
            //fixme implement intersection iterator
            auto barcodes = GetIntersection(first, second);
            size_t result = 0;
            for (auto barcode: barcodes) {
                if (g_.length(first) <= gap_threshold or
                        get_max_pos(first, barcode) > g_.length(first) - gap_threshold) {
                    if (g_.length(second) <= gap_threshold or
                            get_min_pos(second, barcode) < gap_threshold) {
                        ++result;
                    }
                }
            }
            return result;
        }

        //barcode should be present on the edge
        size_t get_min_pos(const EdgeId& edge, int64_t barcode) {
            VERIFY(has_barcode(edge, barcode));
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            auto info_it = entry_it->second.get_barcode(barcode);
            size_t frame_size = entry_it->second.GetFrameSize();
            return info_it->second.GetLeftMost() * frame_size;
        }

        size_t get_max_pos(const EdgeId& edge, int64_t barcode) {
            VERIFY(has_barcode(edge, barcode));
            auto entry_it = mapper_->GetEntryTailsIterator(edge);
            auto info_it = entry_it->second.get_barcode(barcode);
            size_t frame_size = entry_it->second.GetFrameSize();
            return info_it->second.GetRightMost() * frame_size;
        }
    };

} //barcode_index

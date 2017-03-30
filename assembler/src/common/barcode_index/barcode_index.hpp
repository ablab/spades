#pragma once

#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include "io/reads/paired_readers.hpp"
#include <common/assembly_graph/paths/mapping_path.hpp>
#include <common/assembly_graph/core/graph.hpp>
#include "common/pipeline/config_struct.hpp"
#include "common/assembly_graph/index/edge_index_builders.hpp"
#include "common/sequence/range.hpp"

using std::string;
using std::istringstream;
using namespace omnigraph;

namespace barcode_index {
    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    typedef RtSeq Kmer;

    template<class barcode_entry_t>
    class BarcodeIndexBuilder;

    template <class barcode_entry_t>
    class BarcodeIndexInfoExtractor;


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
        std::unordered_map <string, uint64_t> codes_;
    public:
        BarcodeEncoder():
                codes_()
        { }

        void AddBarcode(const string &barcode) {
            auto it = codes_.find(barcode);
            if (it == codes_.end()) {
                size_t encoder_size = codes_.size();
                codes_[barcode] = encoder_size;
            }
        }

        uint64_t GetCode (const string& barcode) const {
            VERIFY(codes_.find(barcode) != codes_.end());
            return codes_.at(barcode);
        }

        size_t GetSize() const {
            return codes_.size();
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

    /**
     * uint64_t wrapper
     */
    class BarcodeId {
        uint64_t int_id_;
        
    public:
        BarcodeId(uint64_t int_id) : int_id_(int_id) {}
        uint64_t int_id() const {
            return int_id_;
        }
        friend bool operator ==(const BarcodeId& left, const BarcodeId& right);
        friend bool operator !=(const BarcodeId& left, const BarcodeId& right);
        friend bool operator <(const BarcodeId& left, const BarcodeId& right);
        friend bool operator >(const BarcodeId& left, const BarcodeId& right);
        friend bool operator <=(const BarcodeId& left, const BarcodeId& right);
        friend bool operator >=(const BarcodeId& left, const BarcodeId& right);
        friend std::ostream& operator<< (std::ostream& stream, const BarcodeId& bid);
    };

    inline bool operator ==(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() == right.int_id();
    }
    inline bool operator !=(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() != right.int_id();
    }
    inline bool operator <(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() < right.int_id();
    }
    inline bool operator >(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() > right.int_id();
    }
    inline bool operator <=(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() <= right.int_id();
    }
    inline bool operator >=(const BarcodeId& left, const BarcodeId& right) {
        return left.int_id() >= right.int_id();
    }
    inline std::ostream& operator<< (std::ostream& stream, const BarcodeId& bid) {
        stream << bid.int_id();
        return stream;
    }

    /**
     This class provides partial interface to BarcodeIndex.
    */
    class AbstractBarcodeIndex {
    public:
    protected:
        const Graph& g_;
    public:
        AbstractBarcodeIndex (const Graph &g) :
                g_(g) {}
        virtual ~AbstractBarcodeIndex() {}

        //Number of entries in the barcode map. Currently equals to the number of edges.
        virtual size_t size() const = 0;

        //Number of barcodes on the beginning/end of the edge
        virtual size_t GetBarcodeNumber(const EdgeId &edge) const = 0;

        //fixme this should be moved to DataScanner
        virtual void ReadEntry(ifstream& fin, const EdgeId& edge) = 0;

        virtual void WriteEntry(ofstream& fin, const EdgeId& edge) = 0;

        //Remove low abundant barcodes
        virtual void Filter(size_t abundancy_threshold, size_t gap_threshold) = 0;

        //Serialize barcode abundancies. Format:
        //abundancy: number of edges.
        virtual bool IsEmpty() = 0;

    };

    /**
     * BarcodeIndex stores information provided by alignment of read clouds to the graph.
     * For every edge we store barcoded reads which are contained on the edge along with additional info.
     * Read cloud is represented by it's barcode
     * The edge contains the cloud if there is a read barcoded by cloud's barcode which is aligned to the edge.
     * Info example: FrameBarcodeInfo
     */
    template <class EdgeEntryT>
    class BarcodeIndex : public AbstractBarcodeIndex {
    friend class BarcodeIndexBuilder<EdgeEntryT>;
    friend class BarcodeIndexInfoExtractor<EdgeEntryT>;
    friend class BarcodeStatisticsCollector;
    protected:
        typedef std::unordered_map <EdgeId, EdgeEntryT> barcode_map_t;
        using AbstractBarcodeIndex::g_;
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
                EdgeEntryT set(*it);
                edge_to_entry_.insert({*it, set});
            }
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

        size_t GetBarcodeNumber(const EdgeId &edge) const override {
            return GetEntry(edge).Size();
        }

        bool IsEmpty() override {
            return size() == 0;
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
            GetEntry(edge).Serialize(fout);
        }

        typename barcode_map_t::const_iterator GetEntryTailsIterator(const EdgeId& edge) const {
            return edge_to_entry_.find(g_.conjugate(edge));
        }

        typename barcode_map_t::const_iterator GetEntryHeadsIterator(const EdgeId& edge) const {
            return edge_to_entry_.find(edge);
        }

        const EdgeEntryT& GetEntry(const EdgeId &edge) const {
            return edge_to_entry_.at(edge);
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

    /**
     * FrameBarcodeInfo approximates the read cloud defined by the barcode and the edge.
     * The edge is split into several bins.
     * Bin is barcoded iff there is at least one barcoded read which aligns to the bin.
     *
     * We store the set of barcoded bins and the number of reads aligned to the edge.
     */
    class FrameBarcodeInfo {
        /**
         * Number of reads aligned to the edge
         */
        size_t count_;
        /**
         * `is_on_frame[i]` is true iff ith bin is barcoded
         */
        boost::dynamic_bitset<> is_on_frame_;
        /**
         * Leftmost barcoded bin
         */
        size_t leftmost_index_;
        /**
         * Rightmost barcoded bin
         */
        size_t rightmost_index_;
    public:

        /**
         *
         * @param frames Number of bin in the edge
         * @return empty info
         */
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
            rightmost_index_ = std::max(rightmost_index_, other.rightmost_index_);
            count_ += other.count_;
        }

        /**
         * @return number of barcoded reads aligned to the edge
         */
        size_t GetCount() const {
            return count_;
        }

        /**
         * @return Leftmost barcoded bin
         */
        size_t GetLeftMost() const {
            return leftmost_index_;
        }

        /**
        * @return Rightmost barcoded bin
        */
        size_t GetRightMost() const {
            return rightmost_index_;
        }

        const boost::dynamic_bitset<>& GetBitSet() const {
            return is_on_frame_;
        }

        /**
         * @param frame index of bin
         * @return true if bin is barcoded, false otherwise
         */
        bool GetFrame(size_t frame) const {
            return is_on_frame_[frame];
        }

        /**
         *
         * @return number of frames
         */
        size_t GetSize() const {
            return is_on_frame_.size();
        }

        /**
         *
         * @return number of barcoded bins
         */
        size_t GetCovered() const {
            return is_on_frame_.count();
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



    template <class entry_info_t>
    class EdgeEntry {
    public:
        typedef std::map <BarcodeId, entry_info_t> barcode_distribution_t;
        typedef entry_info_t barcode_info_t;

    protected:
        EdgeId edge_;
        barcode_distribution_t barcode_distribution_;

    public:
        EdgeEntry():
                edge_(), barcode_distribution_() {};
        EdgeEntry(const EdgeId& edge) :
                edge_(edge), barcode_distribution_() {}

        virtual ~EdgeEntry() {}

        const barcode_distribution_t& GetDistribution() const {
            return barcode_distribution_;
        }

        EdgeId GetEdge() const {
            return edge_;
        }

        //fixme move to extractor
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

        size_t Size() const {
            return barcode_distribution_.size();
        }

        virtual void Serialize(ofstream& fout) const {
            SerializeDistribution(fout);
        }

        virtual void Deserialize(ifstream& fin) {
            DeserializeDistribution(fin);
        }

        typename barcode_distribution_t::const_iterator begin() const {
            return barcode_distribution_.begin();
        }

        typename barcode_distribution_t::const_iterator end() const {
            return barcode_distribution_.end();
        }

        typename barcode_distribution_t::const_iterator cbegin() const {
            return barcode_distribution_.cbegin();
        }

        typename barcode_distribution_t::const_iterator cend() const {
            return barcode_distribution_.cend();
        }

        bool has_barcode(const BarcodeId& barcode) const {
            return barcode_distribution_.find(barcode) != barcode_distribution_.end();
        }

        typename barcode_distribution_t::const_iterator get_barcode(const BarcodeId& barcode) const {
            return barcode_distribution_.find(barcode);
        }

    protected:
        void SerializeDistribution(ofstream &fout) const {
            fout << barcode_distribution_.size() << endl;
            for (auto entry : barcode_distribution_) {
                fout << entry.first.int_id() << ' ' << entry.second << endl;
            }
        }

        void DeserializeDistribution(ifstream &fin) {
            size_t distr_size;
            fin >> distr_size;
            for (size_t i = 0; i < distr_size; ++i) {
                uint64_t int_id;
                entry_info_t info;
                fin >> int_id >> info;
                BarcodeId bid(int_id);
                InsertInfo(bid, info);
            }
        }

        virtual void InsertInfo(const BarcodeId& code, const barcode_info_t& info) = 0;
        virtual void InsertBarcode(const BarcodeId& code, const size_t count, const Range& range) = 0;
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
        void InsertInfo(const BarcodeId& barcode, const SimpleBarcodeInfo &info) {
            if (barcode_distribution_.find(barcode) == barcode_distribution_.end()) {
                barcode_distribution_.insert({barcode, info});
            }
            else {
                barcode_distribution_.at(barcode).Update(info);
            }
        }

        void InsertBarcode(const BarcodeId& barcode, const size_t count, const Range& range) {
            if (barcode_distribution_.find(barcode) == barcode_distribution_.end()) {
                SimpleBarcodeInfo info(count, range);
                barcode_distribution_.insert({barcode, info});
            }
            else {
                barcode_distribution_.at(barcode).Update(count, range);
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

        size_t GetNumberOfFrames() const {
            return number_of_frames_;
        }

    protected:
        void InsertInfo(const BarcodeId& barcode, const FrameBarcodeInfo &info) {
            if (barcode_distribution_.find(barcode) == barcode_distribution_.end()) {
                barcode_distribution_.insert({barcode, info});
            }
            else {
                barcode_distribution_.at(barcode).Update(info);
            }
        }

        void InsertBarcode(const BarcodeId& barcode, const size_t count, const Range& range) {
            if (barcode_distribution_.find(barcode) == barcode_distribution_.end()) {
                FrameBarcodeInfo info(number_of_frames_) ;
                barcode_distribution_.insert({barcode, info});
            }
            else {
                size_t left_frame = GetFrameFromPos(range.start_pos);
                size_t right_frame = GetFrameFromPos(range.end_pos);
                DEBUG("Range: " << range);
                DEBUG("Frames: " << left_frame << " " << right_frame);
                barcode_distribution_.at(barcode).Update(count, left_frame, right_frame);
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
} //barcode_index

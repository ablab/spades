//***************************************************************************
//* Copyright (c) 2017-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/index/edge_index_builders.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "io/binary/binary.hpp"
#include "io/reads/paired_readers.hpp"
#include "sequence/range.hpp"

#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cuckoo/cuckoohash_map.hh>

using std::string;
using std::istringstream;
using namespace omnigraph;

namespace barcode_index {
typedef RtSeq Kmer;

class FrameMapperBuilder;
template <class Graph, class BarcodeEntryT>
class BarcodeIndexInfoExtractor;
typedef uint64_t BarcodeId;

/**
 This class provides partial interface to BarcodeIndex.
*/
template <class Graph>
class AbstractBarcodeIndex {
public:
    typedef typename Graph::EdgeId EdgeId;

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

    virtual void ReadEntry(std::ifstream& fin, const EdgeId& edge) = 0;
    virtual void WriteEntry(std::ofstream& fin, const EdgeId& edge) = 0;

    //Remove low abundant barcodes
    virtual void Filter(size_t abundancy_threshold, size_t gap_threshold) = 0;

    virtual bool IsEmpty() = 0;

};

template <class Graph, class EdgeEntryT>
class ConcurrentBarcodeIndexBuffer {
  public:
    typedef typename Graph::EdgeId EdgeId;
    typedef cuckoohash_map<EdgeId, EdgeEntryT> StorageMap;

    ConcurrentBarcodeIndexBuffer(const Graph &g) : g_(g), edge_to_entry_() {}
    virtual ~ConcurrentBarcodeIndexBuffer() {clear();}

    void clear() {
        edge_to_entry_.clear();
    }

    typename StorageMap::locked_table lock_table() {
        return edge_to_entry_.lock_table();
    }

    void InsertEntry(EdgeEntryT &&entry) {
        edge_to_entry_.insert(entry);
    }

    void InsertBarcode(const BarcodeId &barcode, const EdgeId &edge, size_t count, const Range &range) {
        edge_to_entry_.update_fn(edge,
                                 [&](EdgeEntryT &second) {
                                   second.InsertBarcode(barcode, count, range);
                                 });
    }
  protected:
    const Graph &g_;
    StorageMap edge_to_entry_;
};

/**
 * BarcodeIndex stores information provided by alignment of read clouds to the graph.
 * For every edge we store barcoded reads which are contained on the edge along with additional info.
 * Read cloud is represented by its barcode
 * The edge contains the cloud if there is a read barcoded by cloud's barcode which is aligned to the edge.
 * Info example: FrameBarcodeInfo
 */
template <class Graph, class EdgeEntryT>
class BarcodeIndex: public AbstractBarcodeIndex<Graph> {
friend class BarcodeIndexInfoExtractor<Graph, EdgeEntryT>;

public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef typename omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    typedef std::unordered_map <EdgeId, EdgeEntryT> barcode_map_t;

    BarcodeIndex (const Graph &g) :
            AbstractBarcodeIndex<Graph>(g),
            edge_to_entry_(),
            number_of_barcodes_(0)
    {}

    BarcodeIndex (const BarcodeIndex& other) = default;

    virtual ~BarcodeIndex() {}

    size_t size() const {
        return edge_to_entry_.size();
    }

    bool empty() const {
        return size() == 0;
    }

    typename barcode_map_t::iterator begin() noexcept {
        return edge_to_entry_.begin();
    }
    typename barcode_map_t::iterator end() noexcept {
        return edge_to_entry_.end();
    }
    typename barcode_map_t::iterator begin() const noexcept {
        return edge_to_entry_.begin();
    }
    typename barcode_map_t::iterator end() const noexcept {
        return edge_to_entry_.end();
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

    void ReadEntry (std::ifstream& fin, const EdgeId& edge) override {
        DEBUG("Reading entry")
        DEBUG("Edge: " << edge.int_id());
        DEBUG("Length: " << g_.length(edge));
        edge_to_entry_[edge].Deserialize(fin);
    }
    void WriteEntry (std::ofstream& fout, const EdgeId& edge) override {
        fout << g_.int_id(edge) << std::endl;
        GetEntry(edge).Serialize(fout);
    }

    virtual void BinRead(std::istream &str) {
        using io::binary::BinRead;

        edge_to_entry_.clear();
        size_t size;
        BinRead(str, size);
        for (size_t i = 0; i < size; ++i) {
            EdgeId edge_id = BinRead<uint64_t>(str);
            auto entry = BinRead<EdgeEntryT>(str);
            edge_to_entry_.insert({std::move(edge_id), std::move(entry)});
        }
    }
    virtual void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, edge_to_entry_.size());
        for (const auto &edge_and_entry: edge_to_entry_) {
            BinWrite(str, edge_and_entry.first.int_id(), edge_and_entry.second);
        }
    }

    typename barcode_map_t::const_iterator GetEntryTailsIterator(const EdgeId& edge) const {
        return edge_to_entry_.find(g_.conjugate(edge));
    }
    typename barcode_map_t::const_iterator GetEntryHeadsIterator(const EdgeId& edge) const {
        return edge_to_entry_.find(edge);
    }

    void InsertEntry(const EdgeId &edge, const EdgeEntryT &entry) {
        edge_to_entry_.insert({edge, entry});
    }
    const EdgeEntryT& GetEntry(const EdgeId &edge) const {
        return edge_to_entry_.at(edge);
    }

    void MoveAssign(ConcurrentBarcodeIndexBuffer<Graph, EdgeEntryT> &from) {
        this->edge_to_entry_.clear();
        auto locked_table = from.lock_table();
        for (auto& kvpair : locked_table) {
            this->edge_to_entry_[kvpair.first] = std::move(kvpair.second);
        }
    }

    void SetNumberOfBarcodes(size_t number_of_barcodes) {
        number_of_barcodes_ = number_of_barcodes;
    }
    size_t GetNumberOfBarcodes() {
        return number_of_barcodes_;
    }

    const Graph& GetGraph() const {
        return g_;
    }

 protected:
    using AbstractBarcodeIndex<Graph>::g_;
    barcode_map_t edge_to_entry_;
    size_t number_of_barcodes_;

    DECL_LOGGER("BarcodeIndex");
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
    friend std::ostream& operator <<(std::ostream& os, const SimpleBarcodeInfo& info);
    friend std::istream& operator >>(std::istream& is, SimpleBarcodeInfo& info);
};

inline std::ostream& operator <<(std::ostream& os, const SimpleBarcodeInfo& info)
{
    os << info.count_ << " " << info.range_.start_pos << " " << info.range_.end_pos;
    return os;
}

inline std::istream& operator >>(std::istream& os, SimpleBarcodeInfo& info)
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
public:
    typedef boost::dynamic_bitset<> IsOnFrameT;

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
        TRACE(count_);
        TRACE(other.count_);
        TRACE(is_on_frame_.size());
        TRACE(other.is_on_frame_.size());
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

    const IsOnFrameT& GetBitSet() const {
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

    void SetCount(size_t count) {
        count_ = count;
    }

    void SetLeftMost(size_t index) {
        leftmost_index_ = index;
    }

    void SetRightMost(size_t index) {
        rightmost_index_ = index;
    }

    void SetBitSet(const IsOnFrameT &bitset) {
        is_on_frame_ = bitset;
    }

    void BinRead(std::istream &str) {
        using io::binary::BinRead;
        auto count = BinRead<size_t>(str);
        SetCount(count);

        auto set_positions = BinRead<std::vector<size_t>>(str);
        VERIFY_DEV(set_positions.back() < is_on_frame_.size());
        TRACE("Last position: " << set_positions.back());
        TRACE("Bitset size: " << is_on_frame_.size());
        for (const auto &pos: set_positions) {
            TRACE("Position: " << pos);
            is_on_frame_.set(pos, true);
        }

        SetLeftMost(is_on_frame_.find_first());
        size_t rightmost = 0;
        for (size_t i = is_on_frame_.size() - 1; i > 0; --i) {
            if (is_on_frame_.test(i)) {
                rightmost = i;
                break;
            }
        }
        SetRightMost(rightmost);
        TRACE("Leftmost: " << GetLeftMost());
        TRACE("Rightmost: " << GetRightMost());
    }

    void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, GetCount());

        std::vector<size_t> set_positions;
        size_t current_set_pos = GetBitSet().find_first();
        TRACE("Size: " << GetBitSet().size());
        while (current_set_pos != IsOnFrameT::npos) {
            TRACE("Current set position: " << current_set_pos);
            set_positions.push_back(current_set_pos);
            current_set_pos = GetBitSet().find_next(current_set_pos);
        }
        BinWrite(str, set_positions);
    }


    friend std::ostream& operator <<(std::ostream& os, const FrameBarcodeInfo& info);
    friend std::istream& operator >>(std::istream& is, FrameBarcodeInfo& info);

 private:
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

    DECL_LOGGER("FrameBarcodeInfo");
};

inline std::ostream& operator <<(std::ostream& os, const FrameBarcodeInfo& info)
{
    os << info.count_ << " " << info.is_on_frame_;
    return os;
}

inline std::istream& operator >>(std::istream& is, FrameBarcodeInfo& info)
{
    is >> info.count_;
    is >> info.is_on_frame_;
    info.leftmost_index_ = info.is_on_frame_.find_first();
    size_t rightmost = 0;
    for (size_t i = info.is_on_frame_.size() - 1; i > 0; --i) {
        if (info.is_on_frame_.test(i)) {
            rightmost = i;
            break;
        }
    }
    info.rightmost_index_ = rightmost;
    return is;
}

template <class Graph, class EntryInfoT>
class EdgeEntry {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef std::map <BarcodeId, EntryInfoT> barcode_distribution_t;
    typedef EntryInfoT barcode_info_t;

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

    size_t Size() const {
        return barcode_distribution_.size();
    }

    virtual void Serialize(std::ofstream& fout) const {
        SerializeDistribution(fout);
    }

    virtual void Deserialize(std::ifstream& fin) {
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

    virtual void BinRead(std::istream &str) {
        barcode_distribution_ = io::binary::BinRead<barcode_distribution_t>(str);
    }

    virtual void BinWrite(std::ostream &str) const {
        io::binary::BinWrite(str, barcode_distribution_);
    }

protected:
    void SerializeDistribution(std::ofstream &fout) const {
        fout << barcode_distribution_.size() << std::endl;
        for (auto entry : barcode_distribution_) {
            fout << entry.first << ' ' << entry.second << std::endl;
        }
    }

    void DeserializeDistribution(std::ifstream &fin) {
        size_t distr_size;
        fin >> distr_size;
        for (size_t i = 0; i < distr_size; ++i) {
            uint64_t int_id;
            EntryInfoT info;
            fin >> int_id >> info;
            BarcodeId bid(int_id);
            InsertInfo(bid, info);
        }
    }

    void InsertInfo(const BarcodeId &barcode, const barcode_info_t &info) {
        auto barcode_result = barcode_distribution_.find(barcode);
        if (barcode_result == barcode_distribution_.end()) {
            barcode_distribution_.insert({barcode, info});
        }
        else {
            barcode_result->second.Update(info);
        }
    }
    virtual void InsertBarcode(const BarcodeId &code, const size_t count, const Range &range) = 0;

    EdgeId edge_;
    barcode_distribution_t barcode_distribution_;
};

template<class Graph>
class SimpleEdgeEntry : public EdgeEntry<Graph, SimpleBarcodeInfo> {
    friend class BarcodeIndex<Graph, SimpleEdgeEntry>;
    friend class BarcodeIndexInfoExtractor<Graph, SimpleEdgeEntry>;
protected:
    typedef typename Graph::EdgeId EdgeId;
    using EdgeEntry<Graph, SimpleBarcodeInfo>::barcode_distribution_t;
    using EdgeEntry<Graph, SimpleBarcodeInfo>::barcode_distribution_;
    using EdgeEntry<Graph, SimpleBarcodeInfo>::edge_;

public:
    SimpleEdgeEntry():
        EdgeEntry<Graph, SimpleBarcodeInfo>() {}
    SimpleEdgeEntry(const EdgeId& edge) :
        EdgeEntry<Graph, SimpleBarcodeInfo>(edge) {}

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

template<class Graph>
class FrameEdgeEntry : public EdgeEntry<Graph, FrameBarcodeInfo> {
    friend class BarcodeIndex<Graph, FrameEdgeEntry>;
    friend class FrameMapperBuilder;
    friend class BarcodeIndexInfoExtractor<Graph, FrameEdgeEntry>;
    friend class ConcurrentBarcodeIndexBuffer<Graph, FrameEdgeEntry>;
protected:
    typedef typename Graph::EdgeId EdgeId;
    using EdgeEntry<Graph, FrameBarcodeInfo>::barcode_distribution_t;
    using EdgeEntry<Graph, FrameBarcodeInfo>::barcode_distribution_;
    using EdgeEntry<Graph, FrameBarcodeInfo>::edge_;
    size_t edge_length_;
    size_t frame_size_;
    size_t number_of_frames_;

public:
    FrameEdgeEntry():
        EdgeEntry<Graph, FrameBarcodeInfo>(),
        edge_length_(0),
        frame_size_(0),
        number_of_frames_(0) {}
    FrameEdgeEntry(const EdgeId& edge, size_t edge_length, size_t frame_size) :
        EdgeEntry<Graph, FrameBarcodeInfo>(edge),
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

    void BinRead(std::istream &str) override {
        using io::binary::BinRead;
        edge_length_ = BinRead<size_t>(str);
        frame_size_ = BinRead<size_t>(str);
        number_of_frames_ = BinRead<size_t>(str);

        barcode_distribution_.clear();
        size_t size;
        BinRead(str, size);
        for (size_t i = 0; i < size; ++i) {
            BarcodeId barcode;
            FrameBarcodeInfo info(number_of_frames_);
            BinRead(str, barcode, info);
            barcode_distribution_.insert({std::move(barcode), std::move(info)});
        }
    }

    void BinWrite(std::ostream &str) const override {
        using io::binary::BinWrite;
        BinWrite(str, edge_length_);
        BinWrite(str, frame_size_);
        BinWrite(str, number_of_frames_);

        size_t size = barcode_distribution_.size();
        BinWrite(str, size);
        for (const auto &entry : barcode_distribution_) {
            BinWrite(str, entry.first, entry.second);
        }
    }

protected:
    void InsertBarcode(const BarcodeId& barcode, const size_t count, const Range& range) override {
        DEBUG("Inserting barcode");
        if (barcode_distribution_.find(barcode) == barcode_distribution_.end()) {
            FrameBarcodeInfo info(number_of_frames_);
            barcode_distribution_.insert({barcode, info});
        }
        size_t left_frame = GetFrameFromPos(range.start_pos);
        size_t right_frame = GetFrameFromPos(range.end_pos);
        DEBUG("Range: " << range);
        DEBUG("Frames: " << left_frame << " " << right_frame);
        DEBUG("Count: " << count);
        VERIFY_DEV(barcode_distribution_.find(barcode) != barcode_distribution_.end());
        barcode_distribution_.at(barcode).Update(count, left_frame, right_frame);
    }


    bool IsFarFromEdgeHead(size_t gap_threshold, const FrameBarcodeInfo& info) {
        return info.GetLeftMost() > gap_threshold / frame_size_;
    }

    bool IsLowReadCount(size_t trimming_threshold, const FrameBarcodeInfo& info) {
        return info.GetCount() < trimming_threshold;
    }

    void SetFrameSize(size_t frame_size) {
        frame_size_ = frame_size;
    }

private:
    //fixme last frame is larger than the others
    size_t GetFrameFromPos(size_t pos) {
        return pos / frame_size_;
    }

    DECL_LOGGER("FrameEdgeEntry");
};

template<class Graph>
class FrameConcurrentBarcodeIndexBuffer: public ConcurrentBarcodeIndexBuffer<Graph, FrameEdgeEntry<Graph>> {
  public:
    FrameConcurrentBarcodeIndexBuffer(const Graph &g, size_t frame_size):
        ConcurrentBarcodeIndexBuffer<Graph, FrameEdgeEntry<Graph>>(g), frame_size_(frame_size) {
    }

    void InitialFillMap() {
        VERIFY_DEV(frame_size_ != 0);
        VERIFY_DEV(edge_to_entry_.empty());
        for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
//            FrameEdgeEntry<Graph> entry(edge, g_.length(edge), frame_size_);
            edge_to_entry_.insert(edge, FrameEdgeEntry<Graph>(edge, g_.length(edge), frame_size_));
        }
    }

    size_t GetFrameSize() {
        return frame_size_;
    }

  private:
    using ConcurrentBarcodeIndexBuffer<Graph, FrameEdgeEntry<Graph>>::g_;
    using ConcurrentBarcodeIndexBuffer<Graph, FrameEdgeEntry<Graph>>::edge_to_entry_;
    size_t frame_size_;
};

template<class Graph>
class FrameBarcodeIndex: public BarcodeIndex<Graph, FrameEdgeEntry<Graph>> {
    friend class FrameMapperBuilder;
    friend class BarcodeIndexInfoExtractor<Graph, FrameEdgeEntry<Graph>>;
 public:
    using BarcodeIndex<Graph, FrameEdgeEntry<Graph>>::barcode_map_t;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;

    FrameBarcodeIndex(const Graph &g, size_t frame_size):
        BarcodeIndex<Graph, FrameEdgeEntry<Graph>>(g), frame_size_(frame_size) {
    }

    size_t GetFrameSize() const {
        return frame_size_;
    }

    void SetFrameSize(size_t frame_size) {
        VERIFY_DEV(frame_size_ == 0);
        frame_size_ = frame_size;
    }

    void InitialFillMap() {
        VERIFY_DEV(frame_size_ != 0);
        VERIFY_DEV(edge_to_entry_.empty());
        edge_it_helper helper(g_);
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            FrameEdgeEntry<Graph> entry(*it, g_.length(*it), frame_size_);
            this->InsertEntry(*it, entry);
        }
    }

 private:
    using BarcodeIndex<Graph, FrameEdgeEntry<Graph>>::g_;
    using BarcodeIndex<Graph, FrameEdgeEntry<Graph>>::edge_to_entry_;
    size_t frame_size_;
};
} //barcode_index

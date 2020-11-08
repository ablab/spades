//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "io/binary/binary.hpp"

#include <algorithm>
#include <atomic>
#include <deque>
#include <vector>

namespace path_extend {

class BidirectionalPath;

struct Gap {
    struct Trash {
        uint32_t previous;
        uint32_t current;
        void BinWrite(std::ostream &str) const {
            using io::binary::BinWrite;
            BinWrite(str, previous);
            BinWrite(str, current);
        };

        void BinRead(std::istream &str) {
            using io::binary::BinRead;
            BinRead(str, previous);
            BinRead(str, current);
        }
    };

    using GapSeqType = std::unique_ptr<std::string>;

    GapSeqType gap_seq;
    int gap;
    Trash trash;

//True if gap is resolved not by ordinary procedure but by some sort of magic, and should not be changed.
    bool is_final;

    static constexpr int INVALID_GAP = std::numeric_limits<int>::min();
    static constexpr bool IS_FINAL_DEFAULT = true;
    #define DEFAULT_TRASH {0, 0}

    static Gap INVALID() {
        return Gap(INVALID_GAP);
    }

    //gap is in k+1-mers and does not know about "trash" regions
    Gap(int gap_, Gap::Trash trash_, bool is_final_ = IS_FINAL_DEFAULT, GapSeqType gap_seq_ = nullptr)
        : gap_seq(std::move(gap_seq_))
        , gap(gap_)
        , trash(trash_)
        , is_final(is_final_)
    {}

    explicit Gap(int gap_ = 0, bool is_final_ = IS_FINAL_DEFAULT, GapSeqType gap_seq_ = nullptr)
        : Gap(gap_, DEFAULT_TRASH, is_final_, std::move(gap_seq_))
    {}

    Gap(GapSeqType gap_seq_, int gap_)
        : Gap(gap_, DEFAULT_TRASH, IS_FINAL_DEFAULT, std::move(gap_seq_))
    {}

    Gap(std::string gap_seq_, int gap_)
        : Gap(gap_, DEFAULT_TRASH, IS_FINAL_DEFAULT, MakeGapSeq(std::move(gap_seq_)))
    {}

    Gap(const Gap& other)
        : Gap(other.gap, other.trash, other.is_final, other.CopyGapSeq())
    {};

    Gap& operator=(const Gap& other) {
        gap_seq = other.CopyGapSeq();
        gap = other.gap;
        trash = other.trash;
        is_final = other.is_final;
        return *this;
    };

    Gap(Gap&&) = default;
    Gap& operator=(Gap&&) = default;

    Gap Conjugate() const {
        if (!gap_seq)
            return Gap(gap, {trash.current, trash.previous}, is_final);
        auto gap_seq_complement = *gap_seq;
        std::reverse(gap_seq_complement.begin(), gap_seq_complement.end());
        nucl_str_complement(const_cast<char*>(gap_seq_complement.data()), gap_seq_complement.size());
        return Gap(gap, {trash.current, trash.previous}, is_final, MakeGapSeq(std::move(gap_seq_complement)));
    }

    bool operator==(const Gap &that) const {
        return gap == that.gap &&
                trash.previous == that.trash.previous && trash.current == that.trash.current &&
                is_final == that.is_final;
    }

    bool operator!=(const Gap &that) const noexcept {
        return !(*this == that);
    }

    int Overlap(size_t k) const noexcept {
        return int(k) - gap;
    }

    int OverlapAfterTrim(size_t k) const noexcept {
        return Overlap(k) - trash.current - trash.previous;
    }

    bool NoTrash() const noexcept {
        return trash.current == 0 && trash.previous == 0;
    }

    GapSeqType CopyGapSeq() const {
        return (gap_seq ? std::make_unique<std::string>(*gap_seq) : nullptr);
    }

    static GapSeqType MakeGapSeq(std::string && gap_seq_) {
        return (gap_seq_.empty() ? nullptr : std::make_unique<std::string>(std::move(gap_seq_)));
    }

    void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, gap_seq);
        BinWrite(str, gap);
        BinWrite(str, trash);
        BinWrite(str, is_final);
    };

    void BinRead(std::istream &str) {
        using io::binary::BinRead;
        BinRead(str, gap_seq);
        BinRead(str, gap);
        BinRead(str, trash);
        BinRead(str, is_final);
    }

    #undef DEFAULT_TRASH
};

inline std::ostream& operator<<(std::ostream& os, Gap gap) {
    return os << "[" << gap.gap << ", " << gap.trash.previous << ", " << gap.trash.current << "], final: "<< gap.is_final;
}

class PathListener {
public:
    virtual void FrontEdgeAdded(debruijn_graph::EdgeId e, BidirectionalPath *path, const Gap &gap) = 0;
    virtual void BackEdgeAdded(debruijn_graph::EdgeId e, BidirectionalPath *path, const Gap &gap) = 0;
    virtual void FrontEdgeRemoved(debruijn_graph::EdgeId e, BidirectionalPath *path) = 0;
    virtual void BackEdgeRemoved(debruijn_graph::EdgeId e, BidirectionalPath *path) = 0;
    virtual ~PathListener() {}
};

class GappedPath {
protected:
    using EdgeId = debruijn_graph::EdgeId;
    std::deque<EdgeId> edges_;
    std::deque<Gap> gaps_; // gap0 -> e0 -> gap1 -> e1 -> ... -> gapN -> eN; gap0 = 0

public:
    GappedPath() = default;
    GappedPath(const std::vector<EdgeId>& path)
        : edges_(path.begin(), path.end())
        , gaps_(path.size(), Gap())
    {}

    GappedPath(const GappedPath&) = default;
    GappedPath(GappedPath&&) = default;
    GappedPath& operator=(const GappedPath&) = default;
    GappedPath& operator=(GappedPath&&) = default;

    size_t Size() const noexcept {
        return edges_.size();
    }

    bool Empty() const noexcept {
        return edges_.empty();
    }

    EdgeId operator[](size_t index) const noexcept {
        return edges_[index];
    }

    EdgeId At(size_t index) const {
        return edges_.at(index);
    }

    const Gap& GapAt(size_t index) const noexcept {
        return gaps_[index];
    }

    void SetGapAt(size_t index, const Gap &gap) {
        gaps_[index] = gap;
    }

    EdgeId Back() const noexcept {
        return edges_.back();
    }

    EdgeId Front() const noexcept {
        return edges_.front();
    }
    void PushBack(EdgeId e, Gap gap = Gap()) {
        VERIFY(!edges_.empty() || gap == Gap());
        edges_.push_back(e);
        gaps_.push_back(std::move(gap));
    }

    void PushBack(GappedPath path, Gap gap = Gap()) {
        if (path.Empty())
            return;
        PushBack(path.Front(), std::move(gap));
        std::move(path.edges_.begin()+1, path.edges_.end(), std::back_inserter(edges_));
        std::move(path.gaps_.begin()+1, path.gaps_.end(), std::back_inserter(gaps_));
    }

    void PushBack(const std::vector<EdgeId>& path, Gap gap = Gap()) {
        if (path.empty())
            return;
        PushBack(path.front(), std::move(gap));
        std::copy(path.begin()+1, path.end(), std::back_inserter(edges_));
        gaps_.resize(gaps_.size() + path.size() - 1, Gap());
    }

    void PopBack() noexcept {
        edges_.pop_back();
        gaps_.pop_back();
    }

    void PopBack(size_t count) {
        VERIFY(count <= Size());
        for (size_t i = 0; i < count; ++i)
            PopBack();
    }

    void PushFront(EdgeId e, Gap gap) {
        if (!Empty())
            gaps_[0] = std::move(gap);
        edges_.push_front(e);
        gaps_.push_front(Gap());
    }

    void PopFront() noexcept {
        edges_.pop_front();
        gaps_.pop_front();
        if (!Empty())
            gaps_.front() = Gap();
    }

    void Clear() noexcept {
        edges_.clear();
        gaps_.clear();
    }

    int FindFirst(EdgeId e) const noexcept {
        auto pos = std::find(edges_.begin(), edges_.end(), e);
        if (pos == edges_.end())
            return -1;
        return static_cast<int>(pos - edges_.begin());
    }

    int FindLast(EdgeId e) const noexcept {
        auto pos = std::find(edges_.rbegin(), edges_.rend(), e);
        return static_cast<int>(edges_.rend() - pos - 1);
    }

    bool Contains(EdgeId e) const noexcept {
        return FindFirst(e) != -1;
    }

    std::vector<size_t> FindAll(EdgeId e, size_t start = 0) const {
        VERIFY(start < Size());
        std::vector<size_t> result;
        for (size_t i = start; i < Size(); ++i) {
            if (edges_[i] == e)
                result.push_back(i);
        }
        return result;
    }

    //TODO is it ok not to compare gaps here?
    bool CompareFrom(size_t from, const GappedPath& sample) const noexcept {
        if (from + sample.Size() > Size())
            return false;

        for (size_t i = 0; i < sample.Size(); ++i) {
            if ((*this)[from + i] != sample[i])
                return false;
        }
        return true;
    }

    int FindFirst(const GappedPath& path, size_t from = 0) const noexcept {
        if (path.Size() > Size()) {
            return -1;
        }
        for (size_t i = from; i <= Size() - path.Size(); ++i) {
            if (CompareFrom(i, path))
                return static_cast<int>(i);
        }
        return -1;
    }

    //TODO: Why just naive search?
    int FindLast(const GappedPath& path) const noexcept {
        if (path.Size() > Size()) {
            return -1;
        }
        for (int i = (int) (Size() - path.Size()); i >= 0; --i) {
            if (CompareFrom((size_t) i, path))
                return i;
        }
        return -1;
    }

    bool Equal(const GappedPath& path) const {
        return operator==(path);
    }

    bool operator==(const GappedPath& path) const {
        return Size() == path.Size() && CompareFrom(0, path);
    }

    bool operator!=(const GappedPath& path) const {
        return !operator==(path);
    }

    GappedPath SubPath(size_t from, size_t to) const {
        VERIFY(from <= to && to <= Size());
        GappedPath result;
        if (from < to) {
            result.PushBack(edges_[from], Gap());
            for (size_t i = from + 1; i < to; ++i)
                result.PushBack(edges_[i], gaps_[i]);
        }
        return result;
    }

    GappedPath SubPath(size_t from) const {
        return SubPath(from, Size());
    }

    auto begin() const noexcept -> decltype(edges_.begin()) {
        return edges_.begin();
    }

    auto end() const noexcept -> decltype(edges_.end()) {
        return edges_.end();
    }

    void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, edges_.size());
        for (auto x : edges_)
            BinWrite(str, x);

        BinWrite(str, gaps_.size());
        for (const auto& x : gaps_)
            BinWrite(str, x);
    };

    void BinRead(std::istream &str) {
        using io::binary::BinRead;
        edges_.resize(BinRead<size_t>(str));
        for (auto& x : edges_)
            BinRead(str, x);

        gaps_.resize(BinRead<size_t>(str));
        for (auto& x : gaps_)
            BinRead(str, x);
    }
};

class BidirectionalPath : public PathListener, public GappedPath {
    static std::atomic<uint64_t> path_id_;

    const debruijn_graph::Graph& g_;
    BidirectionalPath* conj_path_;
    // Length from beginning of i-th edge to path end: L(e_i + gap_(i+1) + e_(i+1) + ... + gap_N + e_N)
    std::deque<size_t> cumulative_len_;
    std::vector<PathListener *> listeners_;
    const uint64_t id_;  //Unique ID
    float weight_;
    int cycle_overlapping_; // in edges; [ < 0 ] => is not cycled

public:
    BidirectionalPath(const debruijn_graph::Graph& g)
            : g_(g),
              conj_path_(nullptr),
              id_(path_id_++),
              weight_(1.0),
              cycle_overlapping_(-1) {
    }

    BidirectionalPath(const debruijn_graph::Graph& g, GappedPath path) 
            : BidirectionalPath(g)
    {
        GappedPath::PushBack(std::move(path));
        cumulative_len_.resize(Size(), 0);

        for (size_t i = 0; i < Size(); ++i)
            cumulative_len_.front() += g_.length(edges_[i]) + gaps_[i].gap;

        for (size_t i = 1; i < Size(); ++i)
            cumulative_len_[i] = cumulative_len_[i - 1] - g_.length(edges_[i - 1]) - gaps_[i].gap;
    }

    BidirectionalPath(const debruijn_graph::Graph& g, std::vector<EdgeId> path)
            : BidirectionalPath(g, GappedPath(std::move(path)))
    {}

    BidirectionalPath(const debruijn_graph::Graph& g, EdgeId e)
            : BidirectionalPath(g) {
        PushBack(e);
    }

    BidirectionalPath(const BidirectionalPath& path)
            : GappedPath(path),
              g_(path.g_),
              conj_path_(nullptr),
              cumulative_len_(path.cumulative_len_),
              listeners_(),
              id_(path_id_++),
              weight_(path.weight_),
              cycle_overlapping_(path.cycle_overlapping_) {
    }

    const debruijn_graph::Graph &g() const noexcept {
        return g_;
    }

    void Subscribe(PathListener * listener) {
        listeners_.push_back(listener);
    }

    void SetConjPath(BidirectionalPath* path) noexcept {
        conj_path_ = path;
    }

    const BidirectionalPath* GetConjPath() const noexcept {
        return conj_path_;
    }

    BidirectionalPath* GetConjPath() noexcept {
        return conj_path_;
    }

    void SetWeight(float w) noexcept {
        weight_ = w;
    }

    double GetWeight() const noexcept {
        return weight_;
    }

    const debruijn_graph::Graph& graph() const noexcept {
        return g_;
    }

    /// @brief Check whether (suffix == prefix). If so, set cycle overlapping to new_overlapping.
    /// @returns whether cycle overlapping is set to new_overlapping.
    bool SetCycleOverlapping(int new_overlapping) noexcept {
        if (static_cast<int>(Size()) < new_overlapping)
            return false;

        for (int i = 0; i < new_overlapping; ++i) {
            if (edges_[i] != edges_[Size() - new_overlapping + i] || GapAt(i) != GapAt(Size() - new_overlapping + i))
                return false;
        }

        cycle_overlapping_ = new_overlapping;
        if (conj_path_)
            conj_path_->cycle_overlapping_ = new_overlapping;
        return true;
    }

    int GetCycleOverlapping() const noexcept {
        return cycle_overlapping_;
    }

    bool IsCycle() const noexcept {
        return cycle_overlapping_ >= 0;
    }

    size_t Length() const {
        if (Empty()) {
            return 0;
        }
        VERIFY(gaps_[0].gap == 0);
        return cumulative_len_[0];
    }

    int ShiftLength(size_t index) const {
        return gaps_[index].gap + (int) g_.length(At(index));
    }

    // Length from beginning of i-th edge to path end for forward directed path: L(e1 + e2 + ... + eN)
    size_t LengthAt(size_t index) const noexcept {
        return cumulative_len_[index];
    }

    size_t GetId() const noexcept {
        return id_;
    }

    void PushBack(EdgeId e, Gap gap = Gap()) {
        VERIFY(!edges_.empty() || gap == Gap());
        if (IsCycle()) {
            VERIFY(e == edges_[cycle_overlapping_]);
            ++cycle_overlapping_;
        }
        GappedPath::PushBack(e, std::move(gap));
        IncreaseLengths(g_.length(e), gaps_.back().gap);
        NotifyBackEdgeAdded(e, gaps_.back());
    }

    void PushBack(const BidirectionalPath& path, Gap gap = Gap()) {
        if (path.Size() > 0) {
            VERIFY(path.GapAt(0) == Gap());
            PushBack(path.At(0), std::move(gap));
            for (size_t i = 1; i < path.Size(); ++i)
                PushBack(path.At(i), path.GapAt(i));
        }
    }

    void PushBack(const std::vector<EdgeId>& path, Gap gap = Gap()) {
        VERIFY(!path.empty());
        PushBack(path[0], std::move(gap));
        for (size_t i = 1; i < path.size(); ++i)
            PushBack(path[i]);
    }

    void PopBack() {
        if (edges_.empty())
            return;

        EdgeId e = edges_.back();
        DecreaseLengths();
        GappedPath::PopBack();
        NotifyBackEdgeRemoved(e);
        DecreaseCycleOverlapping();
    }

    void PopBack(size_t count) {
        for (size_t i = 0; i < count; ++i) {
            PopBack();
        }
    }

    void Clear() {
        while (!Empty()) {
            PopBack();
        }
    }

    void FrontEdgeAdded(EdgeId, BidirectionalPath*, const Gap&) override {
        //FIXME is it ok to be empty?
    }

    void BackEdgeAdded(EdgeId e, BidirectionalPath*, const Gap& gap) override {
        PushFront(g_.conjugate(e), gap.Conjugate());
    }

    void FrontEdgeRemoved(EdgeId, BidirectionalPath*) override {
    }

    void BackEdgeRemoved(EdgeId, BidirectionalPath *) override {
        PopFront();
    }

    bool Contains(debruijn_graph::VertexId v) const {
        for(auto edge : edges_) {
            if(g_.EdgeEnd(edge) == v || g_.EdgeStart(edge) == v ) {
                return true;
            }
        }
        return false;
    }

    BidirectionalPath SubPath(size_t from, size_t to) const {
        return BidirectionalPath(g_, GappedPath::SubPath(from, to));
    }

    BidirectionalPath SubPath(size_t from) const {
        return SubPath(from, Size());
    }

    double Coverage() const {
        double cov = 0.0;

        for (size_t i = 0; i < Size(); ++i) {
            cov += g_.coverage(edges_[i]) * (double) g_.length(edges_[i]);
        }
        return cov / (double) Length();
    }

    BidirectionalPath Conjugate() const {
        BidirectionalPath result(g_);
        if (Empty()) {
            return result;
        }
        result.PushBack(g_.conjugate(Back()));
        for (int i = ((int) Size()) - 2; i >= 0; --i) {
            result.PushBack(g_.conjugate(edges_[i]), gaps_[i + 1].Conjugate());
        }
        result.cycle_overlapping_ = cycle_overlapping_;
        return result;
    }

    bool IsCircular() {
        return (Size() > 0 &&
                g_.EdgeStart(Front()) == g_.EdgeEnd(Back()));
    }

    void PrintDEBUG() const;
    void PrintINFO() const;
    void Print(std::ostream &os) const;
    std::string str() const;

private:
    std::vector<std::string> PrintLines() const;

    void IncreaseLengths(size_t length, int gap) {
        for (auto iter = cumulative_len_.begin(); iter != cumulative_len_.end(); ++iter) {
            *iter += length + gap;
        }
        cumulative_len_.push_back(length);
    }

    void DecreaseLengths() {
        size_t length = g_.length(edges_.back()) + gaps_.back().gap;

        for (auto iter = cumulative_len_.begin(); iter != cumulative_len_.end(); ++iter) {
            *iter -= length;
        }
        cumulative_len_.pop_back();
    }

    void NotifyFrontEdgeAdded(EdgeId e, const Gap& gap) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->FrontEdgeAdded(e, this, gap);
        }
    }

    void NotifyBackEdgeAdded(EdgeId e, const Gap& gap) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->BackEdgeAdded(e, this, gap);
        }
    }

    void NotifyFrontEdgeRemoved(EdgeId e) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->FrontEdgeRemoved(e, this);
        }
    }

    void NotifyBackEdgeRemoved(EdgeId e) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->BackEdgeRemoved(e, this);
        }
    }

    void PushFront(EdgeId e, const Gap& gap) {
        if (IsCycle()) {
            VERIFY(e == edges_[Size() - cycle_overlapping_ - 1]);
            ++cycle_overlapping_;
        }

        GappedPath::PushFront(e, gap);

        int length = (int) g_.length(e);
        if (cumulative_len_.empty()) {
            cumulative_len_.push_front(length);
        } else {
            cumulative_len_.push_front(cumulative_len_.front() + length + gap.gap);
        }
        NotifyFrontEdgeAdded(e, gap);
    }

    void PopFront() {
        EdgeId e = edges_.front();
        cumulative_len_.pop_front();
        GappedPath::PopFront();

        NotifyFrontEdgeRemoved(e);
        DecreaseCycleOverlapping();
    }

    void DecreaseCycleOverlapping() noexcept {
        cycle_overlapping_ -= IsCycle();
    }

    DECL_LOGGER("BidirectionalPath");
};

inline int SkipOneGap(debruijn_graph::EdgeId end, const BidirectionalPath& path, int gap, int pos, bool forward) {
    size_t len = 0;
    while (pos < (int) path.Size() && pos >= 0 && end != path.At(pos) && (int) len < 2 * gap) {
        len += path.graph().length(path.At(pos));
        forward ? pos++ : pos--;
    }
    if (pos < (int) path.Size() && pos >= 0 && end == path.At(pos)) {
        return pos;
    }
    return -1;
}

inline void SkipGaps(const BidirectionalPath &path1, size_t &cur_pos1, int gap1,
                     const BidirectionalPath &path2, size_t &cur_pos2, int gap2,
                     bool use_gaps, bool forward) {
    if (use_gaps) {
        if (gap1 > 0 && gap2 <= 0) {
            int temp2 = SkipOneGap(path1.At(cur_pos1), path2, gap1, (int) cur_pos2, forward);
            if (temp2 >= 0) {
                cur_pos2 = (size_t) temp2;
            }
        } else if (gap2 > 0 && gap1 <= 0) {
            int temp1 = SkipOneGap(path2.At(cur_pos2), path1, gap2, (int) cur_pos1, forward);
            if (temp1 >= 0) {
                cur_pos1 = (size_t) temp1;
            }
        } else if (gap1 > 0 && gap2 > 0 && gap1 != gap2) {
            DEBUG("not equal gaps in two paths!!!");
        }
    }
}


//Try do ignore multiple loop traversals
inline size_t FirstNotEqualPosition(const BidirectionalPath &path1, size_t pos1,
                                    const BidirectionalPath &path2, size_t pos2,
                                    bool use_gaps) {
    int cur_pos1 = (int) pos1;
    int cur_pos2 = (int) pos2;
    int gap1 = path1.GapAt(cur_pos1).gap;
    int gap2 = path2.GapAt(cur_pos2).gap;
    while (cur_pos1 >= 0 && cur_pos2 >= 0) {
        if (path1.At(cur_pos1) == path2.At(cur_pos2)) {
            cur_pos1--;
            cur_pos2--;
        } else {
            DEBUG("Not Equal at " << cur_pos1 << " and " << cur_pos2);
            return cur_pos1;
        }
        if (cur_pos1 >= 0 && cur_pos2 >= 0) {
            size_t p1 = (size_t) cur_pos1;
            size_t p2 = (size_t) cur_pos2;
            SkipGaps(path1, p1, gap1, path2, p2, gap2, use_gaps, false);
            cur_pos1 = (int) p1;
            cur_pos2 = (int) p2;
            gap1 = path1.GapAt(cur_pos1).gap;
            gap2 = path2.GapAt(cur_pos2).gap;
        }
    }
    DEBUG("Equal!!");
    return -1UL;
}

inline bool EqualBegins(const BidirectionalPath &path1, size_t pos1,
                        const BidirectionalPath &path2, size_t pos2,
                        bool use_gaps) {
    DEBUG("Checking for equal begins");
    return FirstNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}

inline size_t LastNotEqualPosition(const BidirectionalPath &path1, size_t pos1,
                                   const BidirectionalPath &path2, size_t pos2,
                                   bool use_gaps) {
    size_t cur_pos1 = pos1;
    size_t cur_pos2 = pos2;
    while (cur_pos1 < path1.Size() && cur_pos2 < path2.Size()) {
        if (path1.At(cur_pos1) == path2.At(cur_pos2)) {
            cur_pos1++;
            cur_pos2++;
        } else {
            return cur_pos1;
        }
        int gap1 = cur_pos1 < path1.Size() ? path1.GapAt(cur_pos1).gap : 0;
        int gap2 = cur_pos2 < path2.Size() ? path2.GapAt(cur_pos2).gap : 0;
        SkipGaps(path1, cur_pos1, gap1, path2, cur_pos2, gap2, use_gaps, true);
    }
    return -1UL;
}

inline bool EqualEnds(const BidirectionalPath &path1, size_t pos1,
                      const BidirectionalPath &path2, size_t pos2,
                      bool use_gaps) {
    return LastNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}

inline bool EndsWithInterstrandBulge(const BidirectionalPath &path) {
    if (path.Empty())
        return false;

    const debruijn_graph::Graph &g = path.g();
    debruijn_graph::EdgeId e = path.Back();
    debruijn_graph::VertexId v1 = g.EdgeStart(e);
    debruijn_graph::VertexId v2 = g.EdgeEnd(e);

    return v2 == g.conjugate(v1) &&
            e != g.conjugate(e) &&
            g.OutgoingEdgeCount(v1) == 2 &&
            g.CheckUniqueIncomingEdge(v1);
}

using GappedPathStorage = std::vector<GappedPath>;

using TrustedPathsContainer = std::vector<GappedPathStorage>;

}  // path extend


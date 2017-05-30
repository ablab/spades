//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * bidirectional_path.h
 *
 *  Created on: Nov 14, 2011
 *      Author: andrey
 */
#pragma once

#include <atomic>
#include <boost/algorithm/string.hpp>
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/connected_component.hpp"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

namespace path_extend {

class BidirectionalPath;

struct Gap {
    int gap;
    uint32_t trash_previous;
    uint32_t trash_current;

    static const int INVALID_GAP = std::numeric_limits<int>::min();

    static const Gap& INVALID() {
        static Gap gap = Gap(INVALID_GAP);
        return gap;
    }

    //gap is in k+1-mers and does not know about "trash" regions
    explicit Gap(int gap_ = 0, uint32_t trash_previous_ = 0, uint32_t trash_current_ = 0)
     : gap(gap_), trash_previous(trash_previous_), trash_current(trash_current_)
     { }

    Gap conjugate() const {
        return Gap(gap, trash_current, trash_previous);
    }

    bool operator==(const Gap &that) const {
        return gap == that.gap && trash_previous == that.trash_previous && trash_current == that.trash_current;
    }

    bool operator!=(const Gap &that) const {
        return !(*this == that);
    }

    int overlap(size_t k) const {
        return int(k) - gap;
    }

    int overlap_after_trim(size_t k) const {
        return overlap(k) - trash_current - trash_previous;
    }

    bool NoTrash() const {
        return trash_current == 0 && trash_previous == 0;
    }
};

inline std::ostream& operator<<(std::ostream& os, Gap gap) {
    return os << "[" << gap.gap << ", " << gap.trash_previous << ", " << gap.trash_current << "]";
}

class PathListener {
public:
    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath *path, const Gap &gap) = 0;
    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath *path, const Gap &gap) = 0;
    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath *path) = 0;
    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath *path) = 0;
    virtual ~PathListener() {
    }
};

class BidirectionalPath : public PathListener {
    static std::atomic<uint64_t> path_id_;

    const Graph& g_;
    std::deque<EdgeId> data_;
    BidirectionalPath* conj_path_;
    // Length from beginning of i-th edge to path end: L(e_i + gap_(i+1) + e_(i+1) + ... + gap_N + e_N)
    std::deque<size_t> cumulative_len_;
    std::deque<Gap> gap_len_;  // e0 -> gap1 -> e1 -> ... -> gapN -> eN; gap0 = 0
    std::vector<PathListener *> listeners_;
    const uint64_t id_;  //Unique ID
    float weight_;

public:
    BidirectionalPath(const Graph& g)
            : g_(g),
              data_(),
              conj_path_(nullptr),
              cumulative_len_(),
              gap_len_(),
              listeners_(),
              id_(path_id_++),
              weight_(1.0) {
    }

    BidirectionalPath(const Graph& g, const std::vector<EdgeId>& path)
            : BidirectionalPath(g) {
        for (size_t i = 0; i < path.size(); ++i) {
            PushBack(path[i]);
        }
        RecountLengths();
    }

    BidirectionalPath(const Graph& g, EdgeId startingEdge)
            : BidirectionalPath(g) {
        PushBack(startingEdge);
    }

    BidirectionalPath(const BidirectionalPath& path)
            : g_(path.g_),
              data_(path.data_),
              conj_path_(nullptr),
              cumulative_len_(path.cumulative_len_),
              gap_len_(path.gap_len_),
              listeners_(),
              id_(path_id_++),
              weight_(path.weight_) {
    }

    void Subscribe(PathListener * listener) {
        listeners_.push_back(listener);
    }

//    void Unsubscribe(PathListener * listener) {
//        for (auto it = listeners_.begin(); it != listeners_.end(); ++it) {
//            if (*it == listener) {
//                listeners_.erase(it);
//                break;
//            }
//        }
//    }

    void SetConjPath(BidirectionalPath* path) {
        conj_path_ = path;
    }

    const BidirectionalPath* GetConjPath() const {
        return conj_path_;
    }

    BidirectionalPath* GetConjPath() {
        return conj_path_;
    }

    void SetWeight(float w) {
        weight_ = w;
    }

    double GetWeight() const {
        return weight_;
    }

    size_t Size() const {
        return data_.size();
    }

    const Graph& graph() const {
        return g_;
    }

    bool Empty() const {
        return data_.empty();
    }

    size_t Length() const {
        if (gap_len_.size() == 0 || cumulative_len_.size() == 0) {
            return 0;
        }
        return cumulative_len_[0] + gap_len_[0].gap;
    }

    //TODO iterators forward/reverse
    EdgeId operator[](size_t index) const {
        return data_[index];
    }

    EdgeId At(size_t index) const {
        return data_[index];
    }

    size_t ShiftLength(size_t index) const {
        return gap_len_[index].gap + g_.length(At(index));
    }

    // Length from beginning of i-th edge to path end for forward directed path: L(e1 + e2 + ... + eN)
    size_t LengthAt(size_t index) const {
        return cumulative_len_[index];
    }

    Gap GapAt(size_t index) const {
        return gap_len_[index];
    }

    size_t GetId() const {
        return id_;
    }

    EdgeId Back() const {
        return data_.back();
    }

    EdgeId Front() const {
        return data_.front();
    }

    void PushBack(EdgeId e, const Gap& gap = Gap()) {
        VERIFY(!data_.empty() || gap == Gap());
        data_.push_back(e);
        gap_len_.push_back(gap);
        IncreaseLengths(g_.length(e), gap.gap);
        NotifyBackEdgeAdded(e, gap);
    }

    void PushBack(const BidirectionalPath& path, const Gap& gap = Gap()) {
        if (path.Size() > 0) {
            VERIFY(path.GapAt(0) == Gap());
            PushBack(path.At(0), gap);
            for (size_t i = 1; i < path.Size(); ++i) {
                PushBack(path.At(i), path.GapAt(i));
            }
        }
    }

    void PopBack() {
        if (data_.empty()) {
            return;
        }
        EdgeId e = data_.back();
        DecreaseLengths();
        gap_len_.pop_back();
        data_.pop_back();
        NotifyBackEdgeRemoved(e);
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
        PushFront(g_.conjugate(e), gap.conjugate());
    }

    void FrontEdgeRemoved(EdgeId, BidirectionalPath*) override {
    }

    void BackEdgeRemoved(EdgeId, BidirectionalPath *) override {
        PopFront();
    }

    int FindFirst(EdgeId e) const {
        for (size_t i = 0; i < Size(); ++i) {
            if (data_[i] == e) {
                return (int) i;
            }
        }
        return -1;
    }

    int FindLast(EdgeId e) const {
        for (int i = (int) Size() - 1; i >= 0; --i) {
            if (data_[i] == e) {
                return i;
            }
        }
        return -1;
    }

    bool Contains(EdgeId e) const {
        return FindFirst(e) != -1;
    }

    bool Contains(VertexId v) const {
        for(auto edge : data_) {
            if(g_.EdgeEnd(edge) == v || g_.EdgeStart(edge) == v ) {
                return true;
            }
        }
        return false;
    }

    vector<size_t> FindAll(EdgeId e, size_t start = 0) const {
        vector<size_t> result;
        for (size_t i = start; i < Size(); ++i) {
            if (data_[i] == e) {
                result.push_back(i);
            }
        }
        return result;
    }

    //TODO is it ok not to compare gaps here?
    bool CompareFrom(size_t from, const BidirectionalPath& sample) const {
        if (from + sample.Size() > Size()) {
            return false;
        }

        for (size_t i = 0; i < sample.Size(); ++i) {
            if (At(from + i) != sample[i]) {
                return false;
            }
        }
        return true;
    }

    size_t CommonEndSize(const BidirectionalPath& p) const {
        if (p.Size() == 0) {
            return 0;
        }
        std::vector<size_t> begins = FindAll(p.At(0));
        for (size_t i = 0; i < begins.size(); ++i) {
            size_t it1 = begins[i];
            size_t it2 = 0;
            while (it2 < p.Size() and At(it1) == p.At(it2)) {
                it1++;
                it2++;
                if (it1 == Size()) {
                    return it2;
                }
            }
        }
        return 0;
    }

    int FindFirst(const BidirectionalPath& path, size_t from = 0) const {
        if (path.Size() > Size()) {
            return -1;
        }
        for (size_t i = from; i <= Size() - path.Size(); ++i) {
            if (CompareFrom(i, path)) {
                return (int) i;
            }
        }
        return -1;
    }
//TODO: Why just naive search?
    int FindLast(const BidirectionalPath& path) const {
        if (path.Size() > Size()) {
            return -1;
        }
        for (int i = (int) (Size() - path.Size()); i >= 0; --i) {
            if (CompareFrom((size_t) i, path)) {
                return i;
            }
        }
        return -1;
    }

    bool Equal(const BidirectionalPath& path) const {
        return operator==(path);
    }

    bool operator==(const BidirectionalPath& path) const {
        return Size() == path.Size() && CompareFrom(0, path);
    }

    bool operator!=(const BidirectionalPath& path) const {
        return !operator==(path);
    }

    void CheckConjugateEnd(size_t max_repeat_length) {
        size_t prev_size = 0;
        while (prev_size != Size()) {
            prev_size = Size();
            FindConjEdges(max_repeat_length);
        }
    }

//    size_t GetComponent(const debruijn_graph::ConnectedComponentCounter &component_counter) const {
//        std::unordered_map <size_t, size_t> component_sizes;
//        for (size_t i = 0; i < this->Size(); i++) {
//            auto e = this->At(i);
//            size_t comp_id = component_counter.GetComponent(e);
//            if (component_sizes.find(comp_id) == component_sizes.end())
//                component_sizes[comp_id] = 0;
//            component_sizes[comp_id] += g_.length(e);
//        }
//        size_t ans = 0;
//        size_t maxans = 0;
//        for (auto pp: component_sizes) {
//            if (pp.second > maxans) {
//                ans = pp.first;
//                maxans = pp.second;
//            }
//        }
//        return ans;
//    }

    void FindConjEdges(size_t max_repeat_length) {
        for (size_t begin_pos = 0; begin_pos < Size(); ++begin_pos) {
            size_t begin = begin_pos;
            vector<size_t> conj_pos = FindAll(g_.conjugate(At(begin_pos)), begin + 1);
            for (auto end_pos = conj_pos.rbegin(); end_pos != conj_pos.rend(); ++end_pos) {
                VERIFY(*end_pos < Size());
                size_t end = *end_pos;
                if (end <= begin) {
                    continue;
                }
                while (begin < end && At(begin) == g_.conjugate(At(end))) {
                    begin++;
                    end--;
                }
                DEBUG("Found palindromic fragment from " << begin_pos << " to " << *end_pos);
                PrintDEBUG();
                VERIFY(*end_pos < Size());
                size_t tail_size = Size() - *end_pos - 1;
                size_t head_size = begin_pos;
                size_t palindrom_half_size = begin - begin_pos;
                size_t head_len = Length() - LengthAt(begin_pos);
                size_t tail_len = *end_pos < Size() - 1 ? LengthAt(*end_pos + 1) : 0;
//TODO : this is not true in case of gaps inside the palindrom_len;
                size_t palindrom_len = (size_t) max((int) LengthAt(begin_pos) - (int) LengthAt(begin), 0);
                size_t between = (size_t) max(0, (int) LengthAt(begin) - (int) (end < Size() - 1 ? LengthAt(end + 1) : 0));
                DEBUG("tail len " << tail_len << " head len " << head_len << " palindrom_len "<< palindrom_len << " between " << between);
                if (palindrom_len <= max_repeat_length) {
                    if (palindrom_len < head_len && palindrom_len < tail_len) {
                        DEBUG("too big head and end");
                        continue;
                    }
                    if (between > palindrom_len) {
                        DEBUG("too big part between");
                        continue;
                    }
                }
                bool delete_tail = tail_size < head_size;
                if (tail_size == head_size) {
                    delete_tail = tail_len < head_len;
                }
                if (delete_tail) {
                    PopBack(tail_size + palindrom_half_size);
                    DEBUG("Deleting tail  because of palindrom removal");
                    return;
                } else {
                    GetConjPath()->PopBack(head_size + palindrom_half_size);
                    DEBUG("Deleting head because of palindrom removal");
                    return;
                }
            }
        }
    }

    BidirectionalPath SubPath(size_t from, size_t to) const {
        VERIFY(from <= to && to <= Size());
        BidirectionalPath result(g_);
        for (size_t i = from; i < to; ++i) {
            result.PushBack(data_[i], i == from ? Gap() : gap_len_[i]);
        }
        return result;
    }

    BidirectionalPath SubPath(size_t from) const {
        return SubPath(from, Size());
    }

    double Coverage() const {
        double cov = 0.0;

        for (size_t i = 0; i < Size(); ++i) {
            cov += g_.coverage(data_[i]) * (double) g_.length(data_[i]);
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
            result.PushBack(g_.conjugate(data_[i]), gap_len_[i + 1].conjugate());
        }

        return result;
    }

    vector<EdgeId> ToVector() const {
        return vector<EdgeId>(data_.begin(), data_.end());
    }

//    bool CameToInterstrandBulge() const {
//        if (Empty())
//            return false;
//
//        EdgeId lastEdge = Back();
//        VertexId lastVertex = g_.EdgeEnd(lastEdge);
//
//        if (g_.OutgoingEdgeCount(lastVertex) == 2) {
//            vector<EdgeId> bulgeEdges(g_.out_begin(lastVertex), g_.out_end(lastVertex));
//            VertexId nextVertex = g_.EdgeEnd(bulgeEdges[0]);
//
//            if (bulgeEdges[0] == g_.conjugate(bulgeEdges[1]) && nextVertex == g_.EdgeEnd(bulgeEdges[1]) && g_.CheckUniqueOutgoingEdge(nextVertex)
//                    && *(g_.out_begin(nextVertex)) == g_.conjugate(lastEdge)) {
//
//                DEBUG("Came to interstrand bulge " << g_.int_id(lastEdge));
//                return true;
//            }
//        }
//        return false;
//    }

    bool IsInterstrandBulge() const {
        if (Empty())
            return false;

        EdgeId lastEdge = Back();
        VertexId lastVertex = g_.EdgeEnd(lastEdge);
        VertexId prevVertex = g_.EdgeStart(lastEdge);

        if (g_.OutgoingEdgeCount(prevVertex) == 2 && g_.IncomingEdgeCount(lastVertex) == 2 && g_.CheckUniqueOutgoingEdge(lastVertex)
                && g_.CheckUniqueIncomingEdge(prevVertex) && *(g_.in_begin(prevVertex)) == g_.conjugate(*(g_.out_begin(lastVertex)))) {

            vector<EdgeId> bulgeEdges(g_.out_begin(prevVertex), g_.out_end(prevVertex));
            EdgeId bulgeEdge = bulgeEdges[0] == lastEdge ? bulgeEdges[1] : bulgeEdges[0];

            if (bulgeEdge == g_.conjugate(lastEdge)) {
                DEBUG("In interstrand bulge " << g_.int_id(lastEdge));
                return true;
            }
        }
        return false;
    }

    void PrintDEBUG() const {
        for (const auto& s: PrintLines()) {
            DEBUG(s);
        }
    }

    void PrintINFO() const {
        for (const auto& s: PrintLines()) {
            INFO(s);
        }
    }

    void Print(std::ostream &os) const {
        if (Empty()) {
            return;
        }
        os << "Path " << GetId() << "\n";
        os << "Length " << Length() << "\n";
        os << "Weight " << weight_ << "\n";
        os << "#, edge (length), gap info, total length, total length from start" << "\n";
        for (size_t i = 0; i < Size(); ++i) {
            os << i << ", " << g_.str(At(i))
               << ", " << GapAt(i)
               << ", " << LengthAt(i)
               << ", " << ((Length() < LengthAt(i)) ? 0 : Length() - LengthAt(i)) << "\n";
        }
    }

    std::string str() const {
        stringstream ss;
        Print(ss);
        return ss.str();
    }

    auto begin() const -> decltype(data_.begin()) {
        return data_.begin();
    }

    auto end() const -> decltype(data_.end()) {
        return data_.end();
    }

private:

    vector<std::string> PrintLines() const {
        auto as_str = str();
        boost::trim(as_str);
        std::vector<std::string> result;
        boost::split(result, as_str, boost::is_any_of("\n"), boost::token_compress_on);
        return result;
    }

    void RecountLengths() {
        cumulative_len_.clear();
        size_t currentLength = 0;
        for (auto iter = data_.rbegin(); iter != data_.rend(); ++iter) {
            currentLength += g_.length((EdgeId) *iter);
            cumulative_len_.push_front(currentLength);
        }
    }

    void IncreaseLengths(size_t length, int gap) {
        for (auto iter = cumulative_len_.begin(); iter != cumulative_len_.end(); ++iter) {
            *iter += length + gap;
        }
        cumulative_len_.push_back(length);
    }

    void DecreaseLengths() {
        size_t length = g_.length(data_.back()) + gap_len_.back().gap;

        for (auto iter = cumulative_len_.begin(); iter != cumulative_len_.end(); ++iter) {
            *iter -= length;
        }
        cumulative_len_.pop_back();
    }

    void NotifyFrontEdgeAdded(EdgeId e, Gap gap) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->FrontEdgeAdded(e, this, gap);
        }
    }

    void NotifyBackEdgeAdded(EdgeId e, Gap gap) {
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

    void PushFront(EdgeId e, Gap gap) {
        data_.push_front(e);
        if (gap_len_.size() > 0) {
            VERIFY(gap_len_[0] == Gap());
            gap_len_[0]= gap;
        }
        gap_len_.push_front(Gap(0, 0, 0));

        int length = (int) g_.length(e);
        if (cumulative_len_.empty()) {
            cumulative_len_.push_front(length);
        } else {
            cumulative_len_.push_front(cumulative_len_.front() + length + gap.gap);
        }
        NotifyFrontEdgeAdded(e, gap);
    }

    void PopFront() {
        EdgeId e = data_.front();
        data_.pop_front();
        gap_len_.pop_front();
        cumulative_len_.pop_front();
        if (!gap_len_.empty()) {
            gap_len_.front() = Gap();
        }

        NotifyFrontEdgeRemoved(e);
    }

    DECL_LOGGER("BidirectionalPath");
};

inline int SkipOneGap(EdgeId end, const BidirectionalPath& path, int gap, int pos, bool forward) {
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

inline void SkipGaps(const BidirectionalPath& path1, size_t& cur_pos1, int gap1, const BidirectionalPath& path2, size_t& cur_pos2, int gap2, bool use_gaps,
                     bool forward) {
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
inline size_t FirstNotEqualPosition(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
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
inline bool EqualBegins(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
    DEBUG("Checking for equal begins");
    return FirstNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}

inline size_t LastNotEqualPosition(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
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

inline bool EqualEnds(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
    return LastNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}

typedef std::pair<BidirectionalPath*, BidirectionalPath*> PathPair;

inline bool compare_path_pairs(const PathPair& p1, const PathPair& p2) {
    if (p1.first->Length() != p2.first->Length() || p1.first->Size() == 0 || p2.first->Size() == 0) {
        return p1.first->Length() > p2.first->Length();
    }
    const Graph& g = p1.first->graph();
    return g.int_id(p1.first->Front()) < g.int_id(p2.first->Front());
}

class PathComparator {
public:
    bool operator()(const BidirectionalPath& p1, const BidirectionalPath& p2) const {
        return p1.GetId() < p2.GetId();
    }

    bool operator()(const BidirectionalPath* p1, const BidirectionalPath* p2) const {
        return p1->GetId() < p2->GetId();
    }
};

typedef set<BidirectionalPath*, PathComparator> BidirectionalPathSet;

template<class Value>
using BidirectionalPathMap = map<BidirectionalPath*, Value, PathComparator>;

typedef std::multiset <BidirectionalPath *, PathComparator> BidirectionalPathMultiset;

class PathContainer {
public:

    typedef std::vector<PathPair> PathContainerT;

    class Iterator : public PathContainerT::iterator {
    public:
        Iterator(const PathContainerT::iterator& iter)
                : PathContainerT::iterator(iter) {
        }
        BidirectionalPath* get() const {
            return this->operator *().first;
        }
        BidirectionalPath* getConjugate() const {
            return this->operator *().second;
        }
    };

    class ConstIterator : public PathContainerT::const_iterator {
    public:
        ConstIterator(const PathContainerT::const_iterator& iter)
                : PathContainerT::const_iterator(iter) {
        }
        BidirectionalPath* get() const {
            return this->operator *().first;
        }
        BidirectionalPath* getConjugate() const {
            return this->operator *().second;
        }
    };

    PathContainer() {
    }

    PathContainer(const PathContainer&) = delete;
    PathContainer& operator=(const PathContainer&) = delete;

    PathContainer(PathContainer&&) = default;
    PathContainer& operator=(PathContainer&&) = default;

    BidirectionalPath& operator[](size_t index) const {
        return *(data_[index].first);
    }

    BidirectionalPath* Get(size_t index) const {
        return data_[index].first;
    }

    BidirectionalPath* GetConjugate(size_t index) const {
        return data_[index].second;
    }

    void DeleteAllPaths() {
        for (size_t i = 0; i < data_.size(); ++i) {
            delete data_[i].first;
            delete data_[i].second;
        }
        clear();
    }

    ~PathContainer() {
        DeleteAllPaths();
    }

    size_t size() const {
        return data_.size();
    }

    void clear() {
        data_.clear();
    }

    void reserve(size_t size) {
        data_.reserve(size);
    }

    bool AddPair(BidirectionalPath* p, BidirectionalPath* cp) {
        p->SetConjPath(cp);
        cp->SetConjPath(p);
        p->Subscribe(cp);
        cp->Subscribe(p);
        data_.push_back(std::make_pair(p, cp));
        return true;
    }

    void SortByLength(bool desc = true) {
        std::stable_sort(data_.begin(), data_.end(), [=](const PathPair& p1, const PathPair& p2) {
            return compare_path_pairs(p1, p2) == desc;
        });
    }

    Iterator begin() {
        return Iterator(data_.begin());
    }

    Iterator end() {
        return Iterator(data_.end());
    }


    ConstIterator begin() const {
        return ConstIterator(data_.begin());
    }

    ConstIterator end() const {
        return ConstIterator(data_.end());
    }

    Iterator erase(Iterator iter) {
        return Iterator(data_.erase(iter));
    }

    void print() const {
        for (size_t i = 0; i < size(); ++i) {
            Get(i)->PrintDEBUG();
            GetConjugate(i)->PrintDEBUG();
        }
    }

    //FIXME implement as filter
    void FilterEmptyPaths() {
        DEBUG ("try to delete empty paths");
        for (Iterator iter = begin(); iter != end();) {
            if (iter.get()->Size() == 0) {
                // FIXME: This is trash. PathContainer should own paths
                delete iter.get();
                delete iter.getConjugate();
                iter = erase(iter);
            } else {
                ++iter;
            }
        }
        DEBUG("empty paths are removed");
    }

    //FIXME move out of this class
    void FilterInterstandBulges() {
        DEBUG ("Try to delete paths with interstand bulges");
        for (Iterator iter = begin(); iter != end(); ++iter) {
            if (iter.get()->IsInterstrandBulge()) {
                iter.get()->PopBack();
            }
            if (iter.getConjugate()->IsInterstrandBulge()) {
                iter.getConjugate()->PopBack();
            }
        }
        DEBUG("deleted paths with interstand bulges");
    }

private:
    std::vector<PathPair> data_;

protected:
    DECL_LOGGER("BidirectionalPath");

};

inline pair<size_t, size_t> ComparePaths(size_t start_pos1, size_t start_pos2,
                                         const BidirectionalPath& path1, const BidirectionalPath& path2,
                                         size_t max_diff) {
    path1.PrintDEBUG();
    path2.PrintDEBUG();
    if (start_pos1 >= path1.Size() || start_pos2 >= path2.Size()) {
        return make_pair(start_pos1, start_pos2);
    }
    const Graph& g = path1.graph();
    size_t cur_pos = start_pos1;
    size_t last2 = start_pos2;
    size_t last1 = cur_pos;
    cur_pos++;
    //todo review logic
    size_t diff_len = 0;
    while (cur_pos < path1.Size()) {
        if (diff_len > max_diff) {
            return make_pair(last1, last2);
        }
        EdgeId e = path1[cur_pos];
        vector<size_t> posistions_in_path2 = path2.FindAll(e);
        bool found = false;
        for (size_t pos2 = 0; pos2 < posistions_in_path2.size(); ++pos2) {
            if (posistions_in_path2[pos2] > last2) {
                int path2_segment_len = int(path2.LengthAt(last2)) - int(g.length(path2.At(last2))) -
                    int(path2.LengthAt(posistions_in_path2[pos2]));
                int diff = abs(path1.GapAt(cur_pos).gap - path2_segment_len);

                VERIFY(diff >= 0);
                if (size_t(diff) > max_diff) {
                    break;
                }
                last2 = posistions_in_path2[pos2];
                last1 = cur_pos;
                DEBUG("found " << cur_pos);
                found = true;
                break;
            }
        }
        if (!found) {
            diff_len += g.length(e) + path1.GapAt(cur_pos).gap;
            DEBUG("not found " << cur_pos << " now diff len " << diff_len);
        } else {
            diff_len = 0;
        }
        cur_pos++;
    }
    return make_pair(last1, last2);
}

inline void DeletePaths(BidirectionalPathSet& paths) {
    for (auto i = paths.begin(); i != paths.end(); ++i) {
        delete (*i);
    }
}

inline void DeletePaths(vector<BidirectionalPath*>& paths) {
    for (auto i = paths.begin(); i != paths.end(); ++i) {
        delete (*i);
    }
}

inline void DeleteMapWithPaths(map<EdgeId, BidirectionalPath*> m) {
    for (auto i = m.begin(); i != m.end(); ++i){
        delete i->second;
    }
}

inline bool CheckOverlap(const BidirectionalPath &path1,
                         const BidirectionalPath &path2,
                         size_t overlap) {
    if (path1.Size() < overlap || path2.Size() < overlap)
        return false;
    for (size_t i = 0; i < overlap; ++i)
        if (path1.At(path1.Size() - overlap + i) != path2.At(i))
            return false;
    return true;
}

inline size_t OverlapSize(const BidirectionalPath &path1, const BidirectionalPath &path2) {
    if (path1.Size() == 0)
        return 0;

    size_t max_overlap = 0;
    for (size_t end_pos : path2.FindAll(path1.Back())) {
        size_t overlap = end_pos + 1;
        if (overlap > max_overlap && CheckOverlap(path1, path2, overlap)) {
            max_overlap = overlap;
        }
    }
    return max_overlap;
}

}  // path extend


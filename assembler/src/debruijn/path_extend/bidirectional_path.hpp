//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * bidirectional_path.h
 *
 *  Created on: Nov 14, 2011
 *      Author: andrey
 */

#ifndef BIDIRECTIONAL_PATH_H_
#define BIDIRECTIONAL_PATH_H_

#include "../debruijn_graph.hpp"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

namespace path_extend {

class BidirectionalPath;

class PathListener {

public:
    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap = 0) = 0;

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap = 0) = 0;

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) = 0;

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) = 0;

    virtual ~PathListener() {

    }

};

class LoopDetectorData {

public:
    //Edge and its weight
    typedef std::map <EdgeId, double> AltenativeMap;

protected:
    size_t iteration_;

    AltenativeMap alternatives_;

public:
    LoopDetectorData(size_t i): iteration_(i), alternatives_()  {
    }

    LoopDetectorData(): alternatives_(){
    }

    LoopDetectorData(const LoopDetectorData& d) {
        iteration_ = d.iteration_;
        alternatives_.insert(d.alternatives_.begin(), d.alternatives_.end());
    }

    size_t GetIteration() const {
        return iteration_;
    }

    const AltenativeMap& GetAlternatives() const {
        return alternatives_;
    }

    void AddAlternative(EdgeId e, double w = 1) {
        alternatives_.insert(std::make_pair(e, w));
    }

    void Clear() {
        alternatives_.clear();
        iteration_ = 0;
    }

    bool operator==(const LoopDetectorData& d) const {
        if (alternatives_.size() != d.alternatives_.size()) {
            return false;
        }

        auto iter2 = d.alternatives_.begin();
        for (auto iter1 = alternatives_.begin(); iter2 != d.alternatives_.end() && iter1 != alternatives_.end(); ++iter1, ++iter2) {
            if (iter1->first != iter2->first || iter1->second != iter2->second) {
                return false;
            }
        }

        return true;
    }

protected:
    DECL_LOGGER("BidirectionalPath")
};


class LoopDetector: public PathListener {

protected:
    const Graph& g_;

    size_t currentIteration_;

    LoopDetectorData * current_;

    std::multimap <EdgeId, LoopDetectorData* > data_;

    BidirectionalPath * path_;

public:
    LoopDetector(const Graph& g_, BidirectionalPath * p_);

    void Clear();

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap);

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap);

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path);

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path);


    virtual ~LoopDetector();


    void AddAlternative(EdgeId e, double w = 1);

    void SelectEdge(EdgeId e, double weight = 1);

    size_t LoopEdges(size_t skip_identical_edges, size_t min_cycle_appearences) const;

    size_t LoopLength(size_t skip_identical_edges, size_t min_cycle_appearences) const;

    bool PathIsLoop(size_t edges) const;

    size_t LastLoopCount(size_t skip_identical_edges, size_t min_cycle_appearences) const;

    size_t LastLoopCount(size_t edges) const;

    bool LoopBecameStable() const;

    size_t GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, std::pair<size_t, size_t> iterationRange) const;

    size_t GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, std::pair<size_t, size_t> iterationRange, double coeff) const;

    bool IsCycled(size_t loopLimit, size_t& skip_identical_edges) const;

    size_t EdgesToRemove(size_t skip_identical_edges, bool fullRemoval = false) const;

    void RemoveLoop(size_t skip_identical_edges, bool fullRemoval = true);

    bool EdgeInShortLoop(EdgeId e) const;
    bool PrevEdgeInShortLoop() const;

    void Print() const {
        INFO("== Detector data_ ==");
        for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
            INFO("Edge " << g_.length(iter->first));

            const LoopDetectorData::AltenativeMap& alts = iter->second->GetAlternatives();
            for(auto alt = alts.begin(); alt != alts.end(); ++alt) {
                INFO("Edge " << g_.length(alt->first) << ", weight " << alt->second);
            }
        }
    }

protected:
    DECL_LOGGER("BidirectionalPath")

};


class BidirectionalPath: public PathListener {

protected:
    DECL_LOGGER("BidirectionalPath")


public:
    BidirectionalPath(const Graph& g)
            : g_(g),
              data_(),
              cumulativeLength_(),
              gapLength_(),
              totalLength_(0),
              loopDetector_(g_, this),
              listeners_(),
              weight_(1.0) {
        Init();
    }

    BidirectionalPath(const Graph& g, const std::vector<EdgeId>& path)
            : g_(g),
              data_(),
              cumulativeLength_(),
              gapLength_(),
              totalLength_(0),
              loopDetector_(g_, this),
              listeners_(),
              weight_(1.0) {
        Init();
        if (path.size() != 0) {
            for (size_t i = 0; i < path.size(); ++i) {
                PushBack(path[i]);
            }
            prev_ = path[path.size() - 1];
            now_ = path[path.size() - 1];
            RecountLengths();
        }
    }

    BidirectionalPath(const Graph& g_, EdgeId startingEdge)
            : g_(g_),
              data_(),
              cumulativeLength_(),
              gapLength_(),
              totalLength_(0),
              loopDetector_(g_, this),
              listeners_(),
              weight_(1.0) {
        Init();
        PushBack(startingEdge);
        prev_ = data_.back();
        now_ = data_.back();
    }

    BidirectionalPath(const BidirectionalPath& path)
            : g_(path.g_),
              data_(path.data_),
              cumulativeLength_(path.cumulativeLength_),
              gapLength_(path.gapLength_),
              totalLength_(path.totalLength_),
              loopDetector_(g_, this),
              listeners_(),
              id_(path.id_),
              weight_(path.weight_) {
        Init();
        overlap_ = path.overlap_;
        has_overlaped_begin_ = path.has_overlaped_begin_;
        has_overlaped_end_ = path.has_overlaped_end_;
        prev_ = data_.back();
        now_ = data_.back();
    }


protected:
	const Graph& g_;

	EdgeId prev_, now_;

	// Edges: e1 e2 ... eN
	std::deque <EdgeId> data_;

	BidirectionalPath* conj_path;

	// Length from beginning of i-th edge to path end for forward directed path:
	// L(e1 + e2 + ... + eN) ... L(eN)
	std::deque <size_t> cumulativeLength_;

	// e1 - gap2 - e2 - ... - gapN - eN
	std::deque <int> gapLength_;

	// L(e1 + ... + eN)
	size_t totalLength_;

	// Cycle analyzer
	LoopDetector loopDetector_;

	//Path listeners
	std::set <PathListener *> listeners_;
    size_t id_;//Unique ID in PathContainer
	double weight_;

protected:
	void Verify() {
	    if (cumulativeLength_.empty() && totalLength_ != 0) {
	        INFO("---" << totalLength_);
	    }
	    else if (!cumulativeLength_.empty() && cumulativeLength_[0] + gapLength_[0] != totalLength_) {
	        INFO("||| " << totalLength_ << " !=  " << cumulativeLength_[0] << " + " << gapLength_[0]);
	        PrintInfo();
	    }
	}

    void RecountLengths() {
        Verify();
        cumulativeLength_.clear();
        size_t currentLength = 0;
        for(auto iter = data_.rbegin(); iter != data_.rend(); ++iter) {
            currentLength += g_.length((EdgeId)*iter);
            cumulativeLength_.push_front(currentLength);
        }

        totalLength_ = currentLength;
        Verify();
    }

    void IncreaseLengths(size_t length, size_t gap) {
        Verify();
        for(auto iter = cumulativeLength_.begin(); iter != cumulativeLength_.end(); ++iter) {
            *iter += length + gap;
        }

        cumulativeLength_.push_back(length);
        totalLength_ += length + gap;
        Verify();
    }

    void DecreaseLengths() {
        Verify();
        size_t length = g_.length(data_.back()) + gapLength_.back();
        for(auto iter = cumulativeLength_.begin(); iter != cumulativeLength_.end(); ++iter) {
            *iter -= length;
        }
        cumulativeLength_.pop_back();
        totalLength_ -= length;
        Verify();
    }

    void NotifyFrontEdgeAdded(EdgeId e, int gap) {
        for (auto i = listeners_.begin(); i != listeners_.end(); ++i) {
            (*i)->FrontEdgeAdded(e, this, gap);
        }
    }

    void NotifyBackEdgeAdded(EdgeId e, int gap) {
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

    void PushFront(EdgeId e, int gap = 0) {
        Verify();
        data_.push_front(e);
        if (gapLength_.size() > 0) {
            gapLength_[0] += gap;
        }
        gapLength_.push_front(0);
        int length = (int) g_.length(e);
        if (cumulativeLength_.empty()) {
            cumulativeLength_.push_front(length);
        } else {
            cumulativeLength_.push_front(length + gap + cumulativeLength_.front());
        }
        totalLength_ += length + gap;
        now_ = data_.back();
        NotifyFrontEdgeAdded(e, gap);
        Verify();
    }

    void PopFront() {
        EdgeId e = data_.front();
        int cur_gap = gapLength_.front();
        if (gapLength_.size() > 1) {
            cur_gap += GapAt(1);
            gapLength_[1] = 0;
        }
        data_.pop_front();
        gapLength_.pop_front();

        cumulativeLength_.pop_front();
        totalLength_ -= (g_.length(e) + cur_gap);
        if (data_.size() > 0) {
			now_ = data_.back();
		}else{
			now_ = prev_;
		}
        NotifyFrontEdgeRemoved(e);
        Verify();
    }

    void SafePopFront() {
        PopFront();
    }

public:
    void Subscribe(PathListener * listener) {
        listeners_.insert(listener);
    }

    void Unsubscribe(PathListener * listener) {
        listeners_.erase(listener);
    }

	void SetConjPath(BidirectionalPath* path) {
		conj_path = path;
	}

	BidirectionalPath* GetConjPath() {
		return conj_path;
	}

    const BidirectionalPath* GetConstConjPath() const {
        return conj_path;
    }

    void SetWeight(double w){
    	weight_ = w;
    }

    double GetWeight() const{
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
	    return totalLength_;
	}

	//Access methods
	EdgeId operator[](size_t index) const {
	    return data_[index];
	}

	EdgeId At(size_t index) const {
	    if (index >= data_.size()) {
	        INFO("no data for position " << index << " size " << data_.size());
	    }
	    return data_[index];
	}

	EdgeId ReverseAt(size_t index) const {
        return data_[data_.size() - index - 1];
    }

	size_t LengthAt(size_t index) const {
	    if(index >= cumulativeLength_.size()) {
	        WARN("no length for position " << index <<" size " << cumulativeLength_.size()
	                <<" path size " << cumulativeLength_.size());
	        print_stacktrace();
	        VERIFY(false);
	    }

	    return cumulativeLength_[index];
	}

	int GapAt(size_t index) const {
	    if (index >= gapLength_.size()) {
	        INFO("no gap for position " << index << " size " << gapLength_.size());
	    }
	    return gapLength_[index];
	}

	size_t GetId() const {
	    return id_;
	}

	EdgeId Head() const {
	    return data_.back();
	}

	EdgeId Back() const {
	    return data_.back();
	}

	EdgeId Front() const {
	    return data_.front();
	}

	LoopDetector& getLoopDetector() {
	    return loopDetector_;
	}

	void PushBack(EdgeId e, int gap = 0) {
	    data_.push_back(e);
	    gapLength_.push_back(gap);
	    IncreaseLengths(g_.length(e), gap);
        now_ = e;
	    NotifyBackEdgeAdded(e, gap);
	}


	void PopBack() {
        if (data_.empty()) {
            return;
        }
        EdgeId e = data_.back();
        DecreaseLengths();
	    gapLength_.pop_back();
        data_.pop_back();
        if (data_.size() > 0) {
        	now_ = data_.back();
        }else{
        	now_ = prev_;
        }
	    NotifyBackEdgeRemoved(e);
	}

	void PopBack(size_t count) {
		for (size_t i = 0; i < count; ++i) {
			PopBack();
		}
	}

	void PushBack(const BidirectionalPath& path) {
	    for (size_t i = 0; i < path.Size(); ++i) {
	        PushBack(path.At(i));
	    }
	}

	void SafePopBack() {
	    PopBack();
	}


	void Clear() {
	    while (!Empty()) {
	        PopBack();
	    }
	}

	void SetId(size_t uid) {
	    id_ = uid;
	}

    virtual void FrontEdgeAdded(EdgeId /*e*/, BidirectionalPath * /*path*/, int /*gap*/) {
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * /*path*/, int gap) {
        PushFront(g_.conjugate(e), gap);
    }

    virtual void FrontEdgeRemoved(EdgeId /*e*/, BidirectionalPath * /*path*/) {
    }

    virtual void BackEdgeRemoved(EdgeId /*e*/, BidirectionalPath * /*path*/) {
        PopFront();
    }

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

	size_t CommonEndSize(const BidirectionalPath& sample) const{
		std::vector<size_t> begins;
		for (size_t i = 0; i < Size(); ++i) {
			if (At(i) == sample.At(0)) {
				begins.push_back(i);
			}
		}
		for (size_t i = 0; i < begins.size(); ++i) {
			size_t it1 = begins[i];
			size_t it2 = 0;
			while (it2 < sample.Size() and At(it1) == sample.At(it2)) {
				it1++;
				it2++;
				if (it1 == Size()) {
					return it2;
				}
			}
		}
		return 0;
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

    vector<size_t> FindAll(EdgeId e, size_t start = 0) const {
        vector<size_t> result;
        for (size_t i = start; i < Size(); ++i) {
            if (data_[i] == e) {
                result.push_back(i);
            }
        }
        return result;
    }


	size_t OverlapEndSize(const BidirectionalPath* path2) const {
		if (Size() == 0) {
			return 0;
		}
		int last1 = (int) Size() - 1;
		int max_over = 0;
		vector<size_t> begins2 = path2->FindAll(At(last1));
		for (size_t begin_i = 0; begin_i < begins2.size(); ++begin_i) {
			int begin2 = (int) begins2[begin_i];
			int cur1 = last1;
			while (begin2 > 0 && cur1 > 0
					&& path2->At(begin2 - 1) == At(cur1 - 1)) {
				cur1--;
				begin2--;
			}
			int over = last1 - cur1 + 1;
			if (begin2 == 0 && cur1 > 0 && over > max_over) {
				max_over = over;
			}
		}
		return max_over;
	}

    int FindFirst(const BidirectionalPath& path) const {
        if (path.Size() > Size()) {
            return -1;
        }
        for (size_t i = 0; i <= Size() - path.Size(); ++i) {
            if (CompareFrom(i, path)) {
                return (int) i;
            }
        }
        return -1;
    }

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


    bool Contains(const BidirectionalPath& path) const {
        return FindFirst(path) != -1;
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

    void CheckConjugateEnd() {
        size_t prev_size = 0;
        while (prev_size != Size()) {
            prev_size = Size();
            FindConjEdges();
        }
    }

    void FindConjEdges() {
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

                if (begin >= end) {
                    DEBUG("Found palindromic fragment from " << begin_pos << " to " << *end_pos);
                    Print();
                    VERIFY(*end_pos < Size());
                    size_t tail_size = Size() - *end_pos - 1;
                    size_t head_size = begin_pos;
                    size_t palindrom_half_size = begin - begin_pos;
                    size_t head_len = Length() - LengthAt(begin_pos);
                    size_t tail_len =  *end_pos < Size() - 1 ? LengthAt(*end_pos + 1) : 0;
                    size_t palindrom_len = Length() - head_len - tail_len;
                    if (palindrom_len < head_len && palindrom_len < tail_len) {
                        DEBUG("too big head and end");
                        return;
                    }
                    bool delete_tail = tail_size < head_size;
                    if (tail_size == head_size) {
                        delete_tail = tail_len < head_len;
                    }
                    if (delete_tail) {
                        PopBack(tail_size + palindrom_half_size);
                        FindConjEdges();
                        return;
                    } else {
                        GetConjPath()->PopBack(head_size + palindrom_half_size);
                        FindConjEdges();
                        return;
                    }
                }
            }
        }
    }

    void CheckGrow()
    {
    	prev_ = data_.back();
    	now_ = data_.back();
    }

    bool CheckPrevious()
	{
	    return prev_ != now_;
	}

    BidirectionalPath SubPath(size_t from, size_t to) const {
        BidirectionalPath result(g_);
        for (size_t i = from; i < to; ++i) {
            result.PushBack(data_[i], gapLength_[i]);
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

    BidirectionalPath Conjugate(size_t id = 0) const {
        BidirectionalPath result(g_);
        if (id == 0) {
            result.SetId(id_ % 2 == 0 ? id_ + 1 : id_ - 1);
        } else {
            result.SetId(id);
        }
        if (Empty()) {
            return result;
        }
        result.PushBack(g_.conjugate(Back()), 0);
        for (int i = ((int) Size()) - 2; i >= 0; --i) {
            result.PushBack(g_.conjugate(data_[i]), gapLength_[i + 1]);
        }

        return result;
    }

    vector<EdgeId> ToVector() const {
        return vector<EdgeId>(data_.begin(), data_.end());
    }

    bool IsInterstrandBulge() const {
        EdgeId lastEdge = Back();
        VertexId lastVertex = g_.EdgeEnd(lastEdge);
        VertexId prevVertex = g_.EdgeStart(lastEdge);

        if (g_.OutgoingEdgeCount(prevVertex) == 2 && g_.IncomingEdgeCount(lastVertex) == 2 &&
                g_.CheckUniqueOutgoingEdge(lastVertex) && g_.CheckUniqueIncomingEdge(prevVertex) &&
                *(g_.in_begin(prevVertex)) == g_.conjugate(*(g_.out_begin(lastVertex)))) {

            vector<EdgeId> bulgeEdges(g_.out_begin(prevVertex), g_.out_end(prevVertex));
            EdgeId bulgeEdge = bulgeEdges[0] == lastEdge ? bulgeEdges[1] : bulgeEdges[0];

            if (bulgeEdge == g_.conjugate(lastEdge)) {
                DEBUG("In interstrand bulge " << g_.int_id(lastEdge));
                return true;
            }
        }
        return false;
    }

    void Print() const {
        DEBUG("Path " << id_);
        DEBUG("Length " << totalLength_);
        DEBUG("Weight " << weight_);
        DEBUG("#, edge, length, gap length, total length");

        for(size_t i = 0; i < Size(); ++i) {
        	DEBUG(i << ", " << g_.int_id(At(i)) << ", " << g_.length(At(i)) << ", " << GapAt(i) << ", " << LengthAt(i));
        }
    }

    void PrintInfo() const {
        INFO("Path " << id_);
        INFO("Length " << totalLength_);
        INFO("Weight " << weight_);
        INFO("#, edge, length, gap length, total length");
        for(size_t i = 0; i < Size(); ++i) {
            INFO(i << ", " << g_.int_id(At(i)) << ", " << g_.length(At(i)) << ", " << GapAt(i) << ", " << LengthAt(i));
        }
    }

    void Print(std::ostream& os) {
        if (Empty()) {
            return;
        }
        os << "Path " << GetId() << endl;
        os << "Length " << Length() << endl;
        os << "#, edge, length, total length" << endl;
        for(size_t i = 0; i < Size(); ++i) {
            os << i << ", " << g_.int_id(At(i)) << ", " << g_.length(At(i)) << ", " << LengthAt(i) << endl;
        }
    }

    void SetOverlapedBeginTo(BidirectionalPath* to) {
        if (has_overlaped_begin_) {
            to->SetOverlapBegin();
        }
        SetOverlapBegin();
        to->SetOverlapEnd();
    }

    void SetOverlapedEndTo(BidirectionalPath* to) {
        if (has_overlaped_end_) {
            to->SetOverlapEnd();
        }
        SetOverlapEnd();
        to->SetOverlapBegin();
    }

    void SetOverlap(bool overlap = true) {
        overlap_ = overlap;
        conj_path->overlap_ = overlap;
    }

    bool HasOverlapedBegin() const {
        return has_overlaped_begin_;
    }

    bool HasOverlapedEnd() const {
        return has_overlaped_end_;
    }

    bool IsOverlap() const {
        return overlap_;
    }

private:

    void SetOverlapBegin(bool overlap = true) {
        if (has_overlaped_begin_ != overlap) {
            has_overlaped_begin_ = overlap;
        }
        if (GetConjPath()->has_overlaped_end_ != overlap) {
            GetConjPath()->has_overlaped_end_ = overlap;
        }
    }

    void SetOverlapEnd(bool overlap = true) {
        GetConjPath()->SetOverlapBegin(overlap);
    }

    void Init() {
        Subscribe(&loopDetector_);
        InitOverlapes();
        id_ = 0;
    }

    void InitOverlapes() {
        has_overlaped_begin_ = false;
        has_overlaped_end_ = false;
        overlap_ = false;
    }

    bool has_overlaped_begin_;
    bool has_overlaped_end_;
    bool overlap_;

};

int SkipOneGap(EdgeId end, const BidirectionalPath& path, int gap, int pos,
               bool forward) {
    size_t len = 0;
    while (pos < (int) path.Size() && pos >= 0 && end != path.At(pos)
            && (int) len < 2 * gap) {
        len += path.graph().length(path.At(pos));
        forward ? pos++ : pos--;
    }
    if (pos < (int) path.Size() && pos >= 0 && end == path.At(pos)) {
        return pos;
    }
    return -1;
}

void SkipGaps(const BidirectionalPath& path1, size_t& cur_pos1, int gap1,
              const BidirectionalPath& path2, size_t& cur_pos2, int gap2,
              bool use_gaps, bool forward) {
    if (use_gaps) {
        if (gap1 > 0 && gap2 <= 0) {
            int temp2 = SkipOneGap(path1.At(cur_pos1), path2, gap1,
                                   (int) cur_pos2, forward);
            if (temp2 >= 0) {
                cur_pos2 = (size_t) temp2;
            }
        } else if (gap2 > 0 && gap1 <= 0) {
            int temp1 = SkipOneGap(path2.At(cur_pos2), path1, gap2,
                                   (int) cur_pos1, forward);
            if (temp1 >= 0) {
                cur_pos1 = (size_t) temp1;
            }
        } else if (gap1 > 0 && gap2 > 0 && gap1 != gap2) {
            DEBUG("not equal gaps in two paths!!!");
        }
    }
}

size_t FirstNotEqualPosition(const BidirectionalPath& path1, size_t pos1,
                             const BidirectionalPath& path2, size_t pos2,
                             bool use_gaps) {
    int cur_pos1 = (int) pos1;
    int cur_pos2 = (int) pos2;
    int gap1 = path1.GapAt(cur_pos1);
    int gap2 = path2.GapAt(cur_pos2);
    while (cur_pos1 >= 0 && cur_pos2 >= 0) {
        if (path1.At(cur_pos1) == path2.At(cur_pos2)) {
            cur_pos1--;
            cur_pos2--;
        } else {
            return cur_pos1;
        }
        if (cur_pos1 >= 0 && cur_pos2 >= 0) {
            size_t p1 = (size_t) cur_pos1;
            size_t p2 = (size_t) cur_pos2;
            SkipGaps(path1, p1, gap1, path2, p2, gap2, use_gaps, false);
            cur_pos1 = (int) p1;
            cur_pos2 = (int) p2;
            gap1 = path1.GapAt(cur_pos1);
            gap2 = path2.GapAt(cur_pos2);
        }
    }
    return -1UL;
}
bool EqualBegins(const BidirectionalPath& path1, size_t pos1,
                 const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
    return FirstNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}

size_t LastNotEqualPosition(const BidirectionalPath& path1, size_t pos1,
                            const BidirectionalPath& path2, size_t pos2,
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
        int gap1 = cur_pos1 < path1.Size() ? path1.GapAt(cur_pos1) : 0;
        int gap2 = cur_pos2 < path2.Size() ? path2.GapAt(cur_pos2) : 0;
        SkipGaps(path1, cur_pos1, gap1, path2, cur_pos2, gap2, use_gaps, true);
    }
    return -1UL;
}
bool EqualEnds(const BidirectionalPath& path1, size_t pos1,
               const BidirectionalPath& path2, size_t pos2, bool use_gaps) {
    return LastNotEqualPosition(path1, pos1, path2, pos2, use_gaps) == -1UL;
}
bool PathIdCompare(const BidirectionalPath* p1, const BidirectionalPath* p2) {
    return p1->GetId() < p2->GetId();
}

bool PathCompare(const BidirectionalPath* p1, const BidirectionalPath* p2) {
    if (PathIdCompare(p1, p2) != PathIdCompare(p2, p1))
        return PathIdCompare(p1, p2);
    if (p1->Length() != p2->Length())
        return p1->Length() < p2->Length();
    return p1->Size() < p2->Size();
}

class PathContainer {

public:

    typedef std::pair<BidirectionalPath*, BidirectionalPath*> PathPair;

    typedef std::vector < PathPair > PathContainerT;

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

    PathContainer()
            : path_id_(0) {
    }

    BidirectionalPath& operator[](size_t index) const {
        return *(data_[index].first);
    }

    BidirectionalPath* Get(size_t index) const {
        return data_[index].first;
    }

    BidirectionalPath* GetConjugate(size_t index) const {
        return data_[index].second;
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

    BidirectionalPath* FindConjugate(BidirectionalPath* p) const {
		return p->GetConjPath();
	}

    bool AddPair(BidirectionalPath* p, BidirectionalPath* cp) {
        p->SetConjPath(cp);
        cp->SetConjPath(p);
        p->Subscribe(cp);
        cp->Subscribe(p);
        p->SetId(++path_id_);
        cp->SetId(++path_id_);
        data_.push_back(std::make_pair(p, cp));
        return true;
    }

    void SortByLength() {
        std::stable_sort(data_.begin(), data_.end(), PathPairComparator());
    }

    Iterator begin(){
        return Iterator(data_.begin());
    }

    Iterator end(){
        return Iterator(data_.end());
    }

    Iterator erase(Iterator iter) {
        return Iterator(data_.erase(iter));
    }

    void print() const {
        for (size_t i = 0; i < size(); ++i) {
            Get(i)->Print();
            GetConjugate(i)->Print();
        }
    }

    void CheckSymmetry() const {
        DEBUG("Checking symmetry");
        for (size_t i = 0; i < size(); ++i) {
            if (Get(i)->Conjugate() != *GetConjugate(i)) {
                Get(i)->Print();
                GetConjugate(i)->Print();
            }
        }
    }

    void FilterEmptyPaths(){
    	DEBUG ("try to delete empty paths");
    	for (Iterator iter = begin(); iter != end(); ){
    		if (iter.get()->Size() == 0){
    			iter = erase(iter);
    		} else {
    			++iter;
    		}
    	}
    	DEBUG("empty paths are removed");
    }

    void ResetPathsId() {
        path_id_ = 0;
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i].first->SetId(++path_id_);
            data_[i].second->SetId(++path_id_);
        }
    }

private:

    class PathPairComparator {
    public:

        bool operator()(const PathPair& p1, const PathPair& p2) const {
            bool result = p1.first->Length() > p2.first->Length();
            if (p1.first->Length() == p2.first->Length()) {
                const Graph& g =  p1.first->graph();
                result = g.int_id(p1.first->Front()) < g.int_id(p2.first->Front());
            }
            return result;
        }

        bool operator()(const PathPair* p1, const PathPair* p2) const {
            return operator ()(*p1, *p2);
        }
    };

    std::vector<PathPair> data_;
    size_t path_id_;

protected:
    DECL_LOGGER("BidirectionalPath")

};



LoopDetector::LoopDetector(const Graph& g_, BidirectionalPath * p_): g_(g_), currentIteration_(0), data_(), path_(p_) {
    current_ = new LoopDetectorData(currentIteration_);
}

void LoopDetector::AddAlternative(EdgeId e, double w) {
    current_->AddAlternative(e, w);
}

void LoopDetector::FrontEdgeAdded(EdgeId /*e*/, BidirectionalPath * /*path*/, int /*gap*/) {

}

void LoopDetector::BackEdgeAdded(EdgeId e, BidirectionalPath * /*path*/, int /*gap*/) {
    current_->AddAlternative(e, 1);
    SelectEdge(e);
}

void LoopDetector::FrontEdgeRemoved(EdgeId /*e*/, BidirectionalPath * /*path*/) {

}

void LoopDetector::BackEdgeRemoved(EdgeId e, BidirectionalPath * /*path*/) {
    auto iter = data_.find(e);

    if (iter != data_.end()) {
        iter = data_.upper_bound(e);
        --iter;
        data_.erase(iter);
    }
}

void LoopDetector::SelectEdge(EdgeId e, double /*weight*/) {
    data_.insert(std::make_pair(e, current_));
    current_ = new LoopDetectorData(++currentIteration_);
}


void LoopDetector::Clear() {
    for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
      delete iter->second;
    }

    data_.clear();
    current_->Clear();
}


LoopDetector::~LoopDetector() {
    Clear();
    delete current_;
}


size_t LoopDetector::LoopEdges(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    EdgeId e = path_->Head();

    size_t count = data_.count(e);
    if (count <= 1 || count < min_cycle_appearences * (skip_identical_edges + 1)) {
        return 0;
    }

    auto iter = data_.upper_bound(e);
    --iter;
    size_t loopSize = iter->second->GetIteration();
    for (size_t i = 0; i < skip_identical_edges + 1; ++i) {
        --iter;
    }
    loopSize -= iter->second->GetIteration();

    return loopSize;
}

bool LoopDetector::PathIsLoop(size_t edges) const {
    for (size_t i = 0; i < edges; ++i) {
        EdgeId e = path_->At(i);
        for (int j = (int) path_->Size() - ((int) edges - (int) i); j >= 0; j -= (int) edges) {
            if (path_->operator [](j) != e) {
                return false;
            }
        }
    }
    return true;
}


size_t LoopDetector::LastLoopCount(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    size_t edges = LoopEdges(skip_identical_edges, min_cycle_appearences);
    return LastLoopCount(edges);
}

size_t LoopDetector::LastLoopCount(size_t edges) const {
    if (edges == 0) {
        return 0;
    }

    BidirectionalPath loop = path_->SubPath(path_->Size() - edges);
    size_t count = 0;
    int i = (int) path_->Size() - (int) edges ;
    int delta = - (int) edges;

    while (i >= 0) {
        if (!path_->CompareFrom(i, loop)) {
            break;
        }
        ++count;
        i += delta;
    }

    return count;
}

bool LoopDetector::IsCycled(size_t loopLimit, size_t& skip_identical_edges) const {
    skip_identical_edges = 0;
    size_t loop_count = LastLoopCount(skip_identical_edges, loopLimit);

    while (loop_count > 0) {
        if (loop_count >= loopLimit) {
            return true;
        }
        loop_count = LastLoopCount(++skip_identical_edges, loopLimit);
    }
    return false;
}

size_t LoopDetector::EdgesToRemove(size_t skip_identical_edges,
                                   bool fullRemoval) const {
    size_t edges = LoopEdges(skip_identical_edges, 1);
    size_t count = LastLoopCount(edges);
    bool onlyCycle = PathIsLoop(edges);
    int result;

    if (onlyCycle || path_->Size() <= count * edges) {
        result = (int) path_->Size() - (int) edges;
    } else if (fullRemoval) {
        result = (int) count * (int) edges;
    } else {
        result = (int) (count - 1) * (int) edges;
    }

    return result < 0 ? 0 : result;
}

void LoopDetector::RemoveLoop(size_t skip_identical_edges, bool fullRemoval) {
    size_t toRemove = EdgesToRemove(skip_identical_edges, fullRemoval);
    for(size_t i = 0; i < toRemove; ++i) {
        path_->SafePopBack();
    }
}


bool LoopDetector::LoopBecameStable() const {
    EdgeId e = path_->Head();

    if (data_.count(e) < 2) {
        return false;
    }

    auto iter = data_.upper_bound(e);
    auto last = --iter;
    auto prev = --iter;

    return prev->second == last->second;
}


size_t LoopDetector::GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, std::pair<size_t, size_t> iterationRange) const {
    auto range = data_.equal_range(loopEdge);

    size_t maxIter = 0;
    double maxWeight = 0;
    for (auto iter = range.first; iter != range.second; ++iter) {
        double w = iter->second->GetAlternatives().find(loopExit)->second;
        if (w > maxWeight &&
                iter->second->GetIteration() >= iterationRange.first && iter->second->GetIteration() <= iterationRange.second) {

            maxIter = iter->second->GetIteration();
            maxWeight = w;
        }
    }
    return maxIter;
}

size_t LoopDetector::GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, std::pair<size_t, size_t> iterationRange, double coeff) const {
    auto range = data_.equal_range(loopEdge);

    size_t maxIter = 0;
    for (auto iter = range.first; iter != range.second; ++iter) {
        if (iter->second->GetAlternatives().find(loopExit)->second * coeff > iter->second->GetAlternatives().find(loopEdge)->second && maxIter > iter->second->GetIteration() &&
                iter->second->GetIteration() >= iterationRange.first && iter->second->GetIteration() <= iterationRange.second) {

            maxIter = iter->second->GetIteration();
        }
    }
    return maxIter;
}

bool LoopDetector::EdgeInShortLoop(EdgeId e) const {
	VertexId v = g_.EdgeEnd(e);
	if (g_.OutgoingEdgeCount(v) != 2) {
		return false;
	}
	auto edges = g_.OutgoingEdges(v);
	bool loop_found = false;
	bool exit_found = false;
	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		if (g_.EdgeEnd(*edge) == g_.EdgeStart(e)) {
			loop_found = true;
		} else {
			exit_found = true;
		}
	}
	return loop_found && exit_found;
}

bool LoopDetector::PrevEdgeInShortLoop() const {
    if (path_->Size() <= 1){
    	return false;
    }
	EdgeId e2 = path_->At(path_->Size() - 1);
    EdgeId e1 = path_->At(path_->Size() - 2);
    VertexId v2 = g_.EdgeEnd(e1);
    if (g_.OutgoingEdgeCount(v2) == 2 && g_.EdgeEnd(e2)== g_.EdgeStart(e1) && g_.EdgeEnd(e1)== g_.EdgeStart(e2)) {
        return EdgeInShortLoop(e1);
    }
    return false;
}

} // path extend

#endif /* BIDIRECTIONAL_PATH_H_ */

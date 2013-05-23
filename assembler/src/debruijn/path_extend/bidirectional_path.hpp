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


class CoordPair: public PathListener {

protected:
    const static int DEFAULT = -1;

    int start_;
    int end_;

public:

    CoordPair(): start_(DEFAULT), end_(DEFAULT) {
    }

    CoordPair(int s, int e): start_(s), end_(e) {
    }


    bool In(int c) {
        if (start_ == -1 && end_ == -1) {
            return false;
        }
        return c >= start_ && c <= end_;
    }

    void Set(int s, int e) {
        start_ = s;
        end_ = e;
    }

    void Clear() {
        Set(DEFAULT, DEFAULT);
    }

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        ++start_;
        ++end_;
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {

    }

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        --start_;
        --end_;
        if (start_ < 0) {
            start_ = 0;
        }

        if (end_ < 0) {
            Clear();
        }
    }

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {

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
};


class LoopDetector: public PathListener {

protected:
    Graph& g_;

    size_t currentIteration_;

    LoopDetectorData * current_;

    std::multimap <EdgeId, LoopDetectorData* > data_;

    BidirectionalPath * path_;

public:
    LoopDetector(Graph& g_, BidirectionalPath * p_);

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

    bool EdgeInShortLoop() const;

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

};

class BidirectionalPath: public PathListener {

protected:
	Graph& g_;

	//Unique ID
	size_t id_;

	EdgeId prev_, now_;

	// Edges: e1 e2 ... eN
	std::deque <EdgeId> data_;

	set<BidirectionalPath*> overlaped_begin;
	set<BidirectionalPath*> overlaped_end;
	bool overlap;

	BidirectionalPath* conj_path;

	// Length from beginning of i-th edge to path end for forward directed path:
	// L(e2 + ... + eN) ... L(eN), 0
	// Length from beginning of the path to the end of i-th edge
	// L(e1), L(e1 + e2) ... L(e1 + ... + eN)
	std::deque <size_t> cumulativeLength_;

	// e1 - gap2 - e2 - ... - gapN - eN
	std::deque <int> gapLength_;

	// L(e1 + ... + eN)
	size_t totalLength_;

	// Cycle analyzer
	LoopDetector loopDetector_;

	CoordPair seedCoords_;

	//Path listeners
	std::vector <PathListener *> listeners_;

	double weight;

protected:
	void Verify() {
	    if (cumulativeLength_.empty() && totalLength_ != 0) {
	        INFO("---" << totalLength_);
	    }
	    else if (!cumulativeLength_.empty() && cumulativeLength_[0] != totalLength_) {
	        INFO("|||" << totalLength_);
	    }
	}

    void RecountLengths() {

        cumulativeLength_.clear();
        size_t currentLength = 0;

        for(auto iter = data_.rbegin(); iter != data_.rend(); ++iter) {
            currentLength += g_.length(*iter);
            cumulativeLength_.push_front(currentLength);
        }

        totalLength_ = currentLength;
        Verify();
    }

    void IncreaseLengths(size_t length) {
        Verify();
        for(auto iter = cumulativeLength_.begin(); iter != cumulativeLength_.end(); ++iter) {
            *iter += length;
        }

        cumulativeLength_.push_back(length);
        totalLength_ += length;
        Verify();
    }

    void DecreaseLengths() {
        size_t length = g_.length(data_.back());
        for(auto iter = cumulativeLength_.begin(); iter != cumulativeLength_.end(); ++iter) {
            *iter -= length;
        }

        cumulativeLength_.pop_back();
        totalLength_ -= length;
        Verify();
    }

    void NotifyFrontEdgeAdded(EdgeId e, int gap) {
        for (size_t i = 0; i < listeners_.size(); ++i) {
            listeners_[i]->FrontEdgeAdded(e, this, gap);
        }
    }

    void NotifyBackEdgeAdded(EdgeId e, int gap) {
        for (size_t i = 0; i < listeners_.size(); ++i) {
            listeners_[i]->BackEdgeAdded(e, this, gap);
        }
    }

    void NotifyFrontEdgeRemoved(EdgeId e) {
        for (size_t i = 0; i < listeners_.size(); ++i) {
            listeners_[i]->FrontEdgeRemoved(e, this);
        }
    }

    void NotifyBackEdgeRemoved(EdgeId e) {
        for (size_t i = 0; i < listeners_.size(); ++i) {
            listeners_[i]->BackEdgeRemoved(e, this);
        }
    }

public:
    void PushFront(EdgeId e, int gap = 0) {
        data_.push_front(e);
        gapLength_.push_front(gap);

        int length = g_.length(e);
        if (cumulativeLength_.empty()) {
            cumulativeLength_.push_front(length);
        } else {
            cumulativeLength_.push_front(length + cumulativeLength_.front());
        }
        totalLength_ += length;
        now_ = data_.back();
        NotifyFrontEdgeAdded(e, gap);
        Verify();
    }
protected:
    void PopFront() {
        EdgeId e = data_.front();
        data_.pop_front();
        gapLength_.pop_front();

        cumulativeLength_.pop_front();
        totalLength_ -= g_.length(e);

        if (data_.size() > 0) {
			now_ = data_.back();
		}else{
			now_ = prev_;
		}
        NotifyFrontEdgeRemoved(e);
        Verify();
    }

    void SafePopFront() {
        if (seedCoords_.In(0)) {
            DEBUG("Cannot remove front edge due to seed restrictions");
            return;
        }

        PopFront();
    }

public:
    void Subscribe(PathListener * listener) {
        listeners_.push_back(listener);
    }

    BidirectionalPath(Graph& g): g_(g), data_(), cumulativeLength_(), gapLength_(), totalLength_(0), loopDetector_(g_, this), seedCoords_(), listeners_(){
	    Subscribe(&loopDetector_);
	    Subscribe(&seedCoords_);
	    overlap = false;
	    weight = 1;
	}

	BidirectionalPath(Graph& g, std::vector < EdgeId > path): g_(g), data_(), cumulativeLength_(), gapLength_(), totalLength_(0), loopDetector_(g_, this), seedCoords_(), listeners_() {

		Subscribe(&loopDetector_);
		Subscribe(&seedCoords_);
		overlap = false;
		if (path.size() != 0) {
            id_ = g_.int_id(path[0]);
            for (size_t i = 0; i < path.size(); ++i) {
                Push(path[i]);
            }
            prev_ = path[path.size() - 1];
            now_ = path[path.size() - 1];
            RecountLengths();
		}
		weight = 1;
	}

	BidirectionalPath(const BidirectionalPath& path): g_(path.g_), id_(path.id_), data_(path.data_), cumulativeLength_(path.cumulativeLength_), gapLength_(path.gapLength_), totalLength_(path.totalLength_),
	        loopDetector_(g_, this), seedCoords_(), listeners_() {
	    Subscribe(&loopDetector_);
	    Subscribe(&seedCoords_);
	    overlap = false;
	    prev_ = data_.back();
	    now_ = data_.back();
	    weight = path.getWeight();
	}

	void setConjPath(BidirectionalPath* path) {
		conj_path = path;
	}

	BidirectionalPath* getConjPath() {
		return conj_path;
	}

    void setWeight(double w){
    	weight = w;
    }

    double getWeight() const{
    	return weight;
    }

	void clearOverlapedEnd() {
		overlaped_end.clear();
		conj_path->overlaped_begin.clear();
	}

	void clearOverlapedBegin() {
		overlaped_begin.clear();
		conj_path->overlaped_end.clear();
	}
private:
	void removeOverlapedBegin(BidirectionalPath* p) {
		auto iter = overlaped_begin.find(p);
		if (iter != overlaped_begin.end()) {
			overlaped_begin.erase(iter);
			iter = conj_path->overlaped_end.find(p->conj_path);
			conj_path->overlaped_end.erase(iter);
		}
	}

	void changeOverlapedBegin(BidirectionalPath* from, BidirectionalPath* to) {
		removeOverlapedBegin(from);
		overlaped_begin.insert(to);
		conj_path->overlaped_end.insert(to->conj_path);
	}

	void changeAllOverlapedBegin(BidirectionalPath* to) {
		for (auto iter = overlaped_begin.begin(); iter != overlaped_begin.end(); ++iter) {
			(*iter)->changeOverlapedEnd(this, to);
		}
	}

public:
	void changeOverlapedBeginTo(BidirectionalPath* to) {
		changeAllOverlapedBegin(to);
		to->addOverlapedBegin(getOverlapedBegin());
		clearOverlapedBegin();
		addOverlapedBegin(to);
		to->addOverlapedEnd(this);
	}

	void changeOverlapedEndTo(BidirectionalPath* to) {
		changeAllOverlapedEnd(to);
		to->addOverlapedEnd(getOverlapedEnd());
		clearOverlapedEnd();
		addOverlapedEnd(to);
		to->addOverlapedBegin(this);
	}
private:
	void changeAllOverlapedEnd(BidirectionalPath* to) {
		for (auto iter = overlaped_end.begin(); iter != overlaped_end.end(); ++iter) {
			(*iter)->changeOverlapedBegin(this, to);
		}
	}

	void removeOverlapedEnd(BidirectionalPath* p) {
		auto iter = overlaped_end.find(p);
		if (iter != overlaped_end.end()) {
			overlaped_end.erase(iter);
			iter = conj_path->overlaped_begin.find(p->conj_path);
			conj_path->overlaped_begin.erase(iter);
		}
	}

	void changeOverlapedEnd(BidirectionalPath* from, BidirectionalPath* to) {
		removeOverlapedEnd(from);
		overlaped_begin.insert(to);
		conj_path->overlaped_begin.insert(to->conj_path);
	}
	void addOverlapedBegin(BidirectionalPath* begin) {
		overlaped_begin.insert(begin);
		conj_path->overlaped_end.insert(begin->conj_path);
	}
	void addOverlapedEnd(BidirectionalPath* end) {
		overlaped_end.insert(end);
		conj_path->overlaped_begin.insert(end->conj_path);
	}
public:
	void addOverlapedBegin(set<BidirectionalPath*>& begin) {
		for (auto iter = begin.begin(); iter != begin.end(); ++iter) {
			addOverlapedBegin(*iter);
		}
	}
	void addOverlapedEnd(set<BidirectionalPath*>& end) {
		for (auto iter = end.begin(); iter != end.end(); ++iter) {
			addOverlapedEnd(*iter);
		}
	}

	bool hasOverlapedBegin() {
		return getOverlapedBegin().size() != 0;
	}

	set<BidirectionalPath*>& getOverlapedBegin() {
		return overlaped_begin;
	}

	bool hasOverlapedEnd() {
		return getOverlapedEnd().size() != 0;
	}

	set<BidirectionalPath*>& getOverlapedEnd() {
		return overlaped_end;
	}

	void setOverlap(bool isOverlap_ = true) {
		overlap = isOverlap_;
		if (!conj_path->isOverlap() != isOverlap_) {
			conj_path->setOverlap(true);
		}
	}

	bool isOverlap() const {
		return overlap;
	}


	size_t Size() const {
	    return data_.size();
	}

	Graph& graph() const {
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
	    return data_[index];
	}

	EdgeId ReverseAt(size_t index) const {
        return data_[data_.size() - index - 1];
    }

	size_t LengthAt(size_t index) const {
	    return cumulativeLength_[index];
	}

	int GapAt(size_t index) const {
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

	//Modification methods
	void SetCurrentPathAsSeed() {
	    seedCoords_.Set(0, Size() - 1);
	}

	void PushBack(EdgeId e, int gap = 0) {
	    data_.push_back(e);
	    gapLength_.push_back(gap);
	    IncreaseLengths(g_.length(e));
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

	void SafePopBack() {
        if (seedCoords_.In(Size() - 1)) {
            DEBUG("Cannot remove back edge due to seed restrictions");
            return;
        }

        PopBack();
	}


	void Clear() {
	    seedCoords_.Clear();
	    while (!Empty()) {
	        PopBack();
	    }
	}


	void Push(EdgeId e, int gap = 0) {
        PushBack(e, gap);
	}

	void SetId(size_t uid) {
	    id_ = uid;
	}

    BidirectionalPath(Graph& g_, EdgeId startingEdge): g_(g_), data_(), cumulativeLength_(), gapLength_(), totalLength_(0), loopDetector_(g_, this), seedCoords_(0, 0), listeners_() {
        Subscribe(&loopDetector_);
        Subscribe(&seedCoords_);
        overlap = false;
        Push(startingEdge);
        id_ = g_.int_id(startingEdge);
        prev_ = data_.back();
        now_ = data_.back();
        weight = 1;
    }

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        PushFront(g_.conjugate(e), gap);
    }

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {
    }

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {
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

    int FindFirst(EdgeId e) {
        for (size_t i = 0; i < Size(); ++i) {
            if (data_[i] == e) {
                return i;
            }
        }
        return -1;
    }

    int FindLast(EdgeId e) {
        for (int i = Size(); i > 0; --i) {
            if (data_[i] == e) {
                return i;
            }
        }
        return -1;
    }

    vector<size_t> FindAll(EdgeId e) const{
        vector<size_t> result;
        for (size_t i = 0; i < Size(); ++i) {
            if (data_[i] == e) {
                result.push_back(i);
            }
        }
        return result;
    }

    size_t OverlapEndSize(const BidirectionalPath* path) const {
		if (Size() == 0) {
			return 0;
		}
		int last_index = Size() - 1;
		int max_overlaped_size = 0;
		vector<size_t> begins_in_path = path->FindAll(At(last_index));
		for (size_t begin_index = 0; begin_index < begins_in_path.size(); ++begin_index) {
			int begin_in_path = begins_in_path[begin_index];
			int index_in_current_path = last_index;
			while (begin_in_path > 0 && index_in_current_path > 0
					&& path->At(begin_in_path - 1) == At(index_in_current_path - 1)) {
				index_in_current_path--;
				begin_in_path--;
			}
			int overlaped_size = last_index - index_in_current_path + 1;
			if (begin_in_path == 0) {
				if (index_in_current_path > 0
						&& overlaped_size > max_overlaped_size) {
					max_overlaped_size = overlaped_size;
				}
			}
		}
		return max_overlaped_size;
	}

    bool Contains(const BidirectionalPath& path) const {
        if (path.Size() > Size()) {
            return false;
        }

        for (size_t i = 0; i <= Size() - path.Size(); ++i) {
            if (CompareFrom(i, path)) {
                return true;
            }
        }
        return true;
    }

    bool StartsWith(const BidirectionalPath& path) const {
        return CompareFrom(0, path);
    }

    bool EndsWith(const BidirectionalPath& path) const {
        if (Size() < path.Size()){
            return false;
        }
        return CompareFrom(Size() - path.Size(), path);
    }

    bool operator==(const BidirectionalPath& path) const {
        return Size() == path.Size() && CompareFrom(0, path);
    }

    bool operator!=(const BidirectionalPath& path) const {
        return !operator==(path);
    }

    void CheckConjugateEnd(){
    	size_t begin = 0;
    	size_t end = Size() - 1;
    	while (begin < end && At(begin) == g_.conjugate(At(end))){
    		begin++;
    		end--;
    	}
    	for (size_t i = 0; i < begin ; ++i){
    		PopBack();
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
            result.Push(data_[i], gapLength_[i]);
        }
        return result;
    }

    BidirectionalPath SubPath(size_t from) const {
        return SubPath(from, Size());
    }

    double Coverage() const {
        double cov = 0.0;

        for (size_t i = 0; i < Size(); ++i) {
            cov += g_.coverage(data_[i]) * g_.length(data_[i]);
        }
        return cov / Length();
    }

    BidirectionalPath Conjugate() const {
        BidirectionalPath result(g_);
        for (size_t i = 0; i < Size(); ++i) {
            result.PushFront(g_.conjugate(data_[i]));
        }
        return result;
    }

    vector<EdgeId> ToVector() const {
        return vector<EdgeId>(data_.begin(), data_.end());
    }

    void Print() const {
        DEBUG("Path " << id_);
        DEBUG("Length " << totalLength_);
        DEBUG("#, edge, length, total length");
        for(size_t i = 0; i < Size(); ++i) {
        	DEBUG(i << ", " << g_.int_id(At(i)) << ", " << g_.length(At(i)) << ", " << LengthAt(i));
        }
    }

    void PrintInfo() const {
        INFO("Path " << id_);
        INFO("Length " << totalLength_);
        INFO("#, edge, length, total length");
        for(size_t i = 0; i < Size(); ++i) {
            INFO(i << ", " << g_.int_id(At(i)) << ", " << g_.length(At(i)) << ", " << GapAt(i));
        }
    }

    void Print(std::ostream& os) {

        os << "Path " << GetId() << endl;

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
};






class PathContainer {

public:

    typedef std::pair<BidirectionalPath*, BidirectionalPath*> PathPair;

    typedef std::vector < PathPair > PathContainerT;

public:

    class Iterator: public PathContainerT::iterator {

    public:
        Iterator(const PathContainerT::iterator& iter): PathContainerT::iterator(iter) {
        }

        BidirectionalPath* get() const {
            return this->operator *().first;
        }

        BidirectionalPath* getConjugate() const {
            return this->operator *().second;
        }
    };

private:
	std::vector < PathPair > data_;

	class PathPairComparator {
	public:

	    bool operator() (const PathPair& p1, const PathPair& p2) const {
	        return p1.first->Length() > p2.first->Length();
	    }

	    bool operator() (const PathPair* p1, const PathPair* p2) const {
	        return p1->first->Length() > p2->first->Length();
	    }
	};

public:
	PathContainer() {
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

    BidirectionalPath* FindConjugate(BidirectionalPath* p) {
		return p->getConjPath();
	}

    bool AddPair(BidirectionalPath* p, BidirectionalPath* cp) {
    	p->setConjPath(cp);
    	cp->setConjPath(p);
    	p->Subscribe(cp);
        cp->Subscribe(p);

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

    void SubscribeAll(PathListener * listener) {
        for (auto iter = begin(); iter != end(); ++iter) {
            iter.get()->Subscribe(listener);
            iter.getConjugate()->Subscribe(listener);
        }
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
};



LoopDetector::LoopDetector(Graph& g_, BidirectionalPath * p_): g_(g_), currentIteration_(0), data_(), path_(p_) {
    current_ = new LoopDetectorData(currentIteration_);
}

void LoopDetector::AddAlternative(EdgeId e, double w) {
    current_->AddAlternative(e, w);
}

void LoopDetector::FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {

}

void LoopDetector::BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
    current_->AddAlternative(e, 1);
    SelectEdge(e);
}

void LoopDetector::FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {

}

void LoopDetector::BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {
    auto iter = data_.find(e);

    if (iter != data_.end()) {
        iter = data_.upper_bound(e);
        --iter;
        data_.erase(iter);
    }
}

void LoopDetector::SelectEdge(EdgeId e, double weight) {
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

//size_t LoopDetector::LoopLength(size_t skip_identical_edges, size_t min_cycle_appearences) const {
//    if (skip_identical_edges > 0) {
//        INFO("loop length " << skip_identical_edges);
//    }
//    size_t edges = LoopEdges(skip_identical_edges, min_cycle_appearences);
//    size_t length = 0;
//
//    for (int i = path_->Size() - edges; i < (int) path_->Size(); ++i) {
//        length += g_.length(path_->At(i));
//    }
//
//    return length;
//}


bool LoopDetector::PathIsLoop(size_t edges) const {
    for (size_t i = 1; i <= edges; ++i) {
        EdgeId e = path_->operator [](path_->Size() - i);
        for (int j = path_->Size() - i - edges; j >= 0; j -= edges) {
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
    int i = path_->Size() - edges ;
    int delta = -edges;

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

size_t LoopDetector::EdgesToRemove(size_t skip_identical_edges, bool fullRemoval) const {
    size_t edges = LoopEdges(skip_identical_edges, 1);
    size_t count = LastLoopCount(edges);
    bool onlyCycle = PathIsLoop(edges);
    int result;

    if (onlyCycle || path_->Size() <= count * edges + 1) {
        result = path_->Size() - edges - 1;
    }
    else if (fullRemoval) {
        result = count * edges - 1;
    } else {
        result = (count - 1) * edges - 1;
    }

    return result < 0 ? 0 : result;
}

void LoopDetector::RemoveLoop(size_t skip_identical_edges, bool fullRemoval) {
    auto toRemove = EdgesToRemove(skip_identical_edges, fullRemoval);
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

bool LoopDetector::EdgeInShortLoop() const {
    EdgeId e = path_->Head();
    VertexId v = g_.EdgeEnd(e);

    if (g_.OutgoingEdgeCount(v) != 2) {
        return false;
    }

    auto edges = g_.OutgoingEdges(v);
    for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
        if (g_.EdgeEnd(*edge) == g_.EdgeStart(e)) {
            return true;
        }
    }

    return false;
}

} // path extend

#endif /* BIDIRECTIONAL_PATH_H_ */

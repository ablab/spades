/*
 * bidirectional_path.h
 *
 *  Created on: Nov 14, 2011
 *      Author: andrey
 */

#ifndef BIDIRECTIONAL_PATH_H_
#define BIDIRECTIONAL_PATH_H_

#include "../new_debruijn.hpp"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;
using debruijn_graph::VerexId;

class LoopDetecor;

class BidirectionalPath {

protected:
	Graph& g;

	std::deque<EdgeId> data;
	std::deque<size_t> cumulativeLength;
	std::deque<size_t> gapLength;

	size_t totalLength;
	bool direction;

	LoopDetector loopDetector;

protected:
	void recountLengths();

	void increaseLengths(size_t length, bool direct = direction);

	void decreaseLengths(bool direct = direction);

public:
	BidirectionalPath(Graph g_): g(g_), data(), cumulativeLength(), gapLength(), direction(true), loopDetector(g_, *this) {
	}

	size_t size() const;

	size_t length() const;

	//Access methods
	EdgeId operator[](size_t index) const;

	EdgeId at(size_t index) const;

	size_t lengthAt(size_t index) const;

	size_t gapAt(size_t index) const;

	bool getDirection() const;

	EdgeId head() const;

	EdgeId back() const;

	EdgeId front() const;

	//Modification methods
	void setDirection(bool d);

	void pushBack(EdgeId e, size_t gap = 0);

	void pushFront(EdgeId e, size_t gap = 0);

	void popBack();

	void popFront();

	void push(EdgeId e, size_t gap = 0);

	bool pop();
};


class PathContainer {
private:
	std::vector<BidirectionalPath*> data;

public:
	PathContainer() {
		data.clear();
	}

	BidirectionalPath* get(size_t index) const;

	BidirectionalPath* getConjugate(size_t index) const;

	void reserve(size_t size);

	bool addPair(BidirectionalPath* p, BidirectionalPath* cp);

	bool addPair(BidirectionalPath& p, BidirectionalPath& cp);
};


#endif /* BIDIRECTIONAL_PATH_H_ */

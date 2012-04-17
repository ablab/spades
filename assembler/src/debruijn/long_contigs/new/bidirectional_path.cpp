//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * bidirectional_path.hpp
 *
 *  Created on: Sep 24, 2011
 *      Author: andrey
 */


#include "bidirectional_path.h"

namespace long_contigs {

void BidirectionalPath::recountLengths() {
	cumulativeLength.clear();
	size_t currentLength = 0;

	if (direction) {
		for(auto iter = data.rbegin(); iter != data.rend(); ++iter) {
			currentLength += g.length(*iter);
			cumulativeLength.push_front(currentLength);
		}
	} else {
		for(auto iter = data.begin(); iter != data.end(); ++iter) {
			cumulativeLength.push_back(currentLength);
			currentLength += g.length(*iter);
		}
	}

	totalLength = currentLength;
}

void BidirectionalPath::increaseLengths(size_t length, bool direct) {
	for(auto iter = cumulativeLength.begin(); iter != cumulativeLength.end(); ++iter) {
		*iter += length;
	}

	if (direct) {
		cumulativeLength.push_back(len);
	} else {
		cumulativeLength.push_front(0);
	}

	totalLength += length;
}

void BidirectionalPath::decreaseLengths(bool direct) {
	size_t length = direct ? g.length(data.back()) : g.length(data.front());
	for(auto iter = cumulativeLength.begin(); iter != cumulativeLength.end(); ++iter) {
		*iter -= length;
	}

	if (direct) {
		cumulativeLength.pop_back();
	} else {
		cumulativeLength.pop_front();
	}

	totalLength -= length;
}



size_t BidirectionalPath::size() const {
	return data.size();
}

size_t BidirectionalPath::length() const {
	return totalLength;
}

//Access methods
EdgeId BidirectionalPath::operator[](size_t index) const {
	return data[index];
}

EdgeId BidirectionalPath::at(size_t index) const {
	return data[index];
}

size_t BidirectionalPath::lengthAt(size_t index) const {
	return cumulativeLength[index];
}

size_t BidirectionalPath::gapAt(size_t index) const {
	return gapLength[index];
}

bool BidirectionalPath::getDirection() const {
	return direction;
}

EdgeId BidirectionalPath::head() const {
	return direction ? data.back() : data.front();
}

EdgeId BidirectionalPath::back() const {
	return data.back();
}

EdgeId BidirectionalPath::front() const {
	return data.front();
}

//Modification methods
void BidirectionalPath::setDirection(bool d) {
	direction = d;
}

void BidirectionalPath::pushBack(EdgeId e, size_t gap) {
	data.push_back(e);
	gaps.push_back(gap);
	increaseLengths(g.length(e), true);
}

void BidirectionalPath::pushFront(EdgeId e, size_t gap) {
	data.push_front(e);
	gaps.push_front(gap);
	increaseLengths(g.length(e), false);
}

void BidirectionalPath::popBack() {
	data.pop_back();
	gaps.pop_back();
	decreaseLengths(true);
}

void BidirectionalPath::popFront() {
	data.pop_front();
	gaps.pop_front();
	decreaseLengths(false);
}


void BidirectionalPath::push(EdgeId e, size_t gap) {
	if (direction) {
		pushBack(e, gap);
	} else {
		pushFront(e, gap);
	}
}

bool BidirectionalPath::pop() {
	if (data.empty()) {
		return false;
	}
	if (direction) {
		popBack();
	} else {
		popFront();
	}
	return true;
}


BidirectionalPath* PathContainer::get(size_t index) const {
	return data[index];
}

BidirectionalPath* PathContainer::getConjugate(size_t index) const {
	return data[(index & 1 == 0) ? index + 1 : index - 1];
}

void PathContainer::reserve(size_t size) {
	data.reserve(size);
}

bool PathContainer::addPair(BidirectionalPath* p, BidirectionalPath* cp) {
	if (p->size() != cp->size() || p->length() != cp->length()) {
		return false;
	}

	data.push_back(p);
	data.push_back(cp);
	return true;
}

bool PathContainer::addPair(BidirectionalPath& p, BidirectionalPath& cp) {
	if (p->size() != cp->size() || p->length() != cp->length()) {
		return false;
	}
	BidirectionalPath np = new BidirectionalPath(p);
	BidirectionalPath ncp = new BidirectionalPath(cp);

	data.push_back(np);
	data.push_back(ncp);
	return true;
}


}


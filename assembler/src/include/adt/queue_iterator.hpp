//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef QUEUE_ITERATOR_HPP_
#define QUEUE_ITERATOR_HPP_

#include "verify.hpp"
#include <set>

template<typename Key, typename Comparator>
class erasable_priority_queue {
private:
	std::set<Key, Comparator> storage_;
public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with VERIFY(false) inside it.
	 */
	erasable_priority_queue(const Comparator& comparator = Comparator()) :
		storage_(comparator) {
	}

	template<typename InputIterator>
	erasable_priority_queue(InputIterator begin, InputIterator end,
			const Comparator& comparator = Comparator()) :
		storage_(begin, end, comparator) {
	}

	void pop() {
		VERIFY(!storage_.empty());
		storage_.erase(storage_.begin());
	}

	const Key& top() const {
		VERIFY(!storage_.empty());
		return *(storage_.begin());
	}

	void push(const Key& key) {
		storage_.insert(key);
	}

	bool erase(const Key& key) {
		bool res = storage_.erase(key) > 0;
		return res;
	}

	bool empty() const {
		return storage_.empty();
	}

	size_t size() const {
		return storage_.size();
	}

	template <class InputIterator>
	void insert ( InputIterator first, InputIterator last ) {
		storage_.insert(first, last);
	}

};

template<typename ElementId, typename Comparator = std::less<ElementId> >
class QueueIterator {
private:
	bool current_actual_;
	bool current_deleted_;
	ElementId current_;
	erasable_priority_queue<ElementId, Comparator> queue_;
protected:

	template<typename InputIterator>
	void insert(InputIterator begin, InputIterator end) {
		queue_.insert(begin, end);
	}

	QueueIterator(const Comparator& comparator = Comparator()) :
		current_actual_(false), current_deleted_(false), queue_(comparator) {
	}

	void push(const ElementId& toAdd) {
		queue_.push(toAdd);
	}

public:
	void erase(const ElementId& toRemove) {
		if (current_actual_ && toRemove == current_) {
			current_deleted_ = true;
		}
		queue_.erase(toRemove);
	}

	bool IsEnd() const {
		return queue_.empty();
	}

	ElementId operator*() {
		VERIFY(!queue_.empty());
		if(!current_actual_ || current_deleted_) {
			current_ = queue_.top();
			current_actual_ = true;
			current_deleted_ = false;
		}
		return current_;
	}

	void operator++() {
		if (!current_actual_) {
			queue_.pop();
		} else if (!current_deleted_) {
			queue_.erase(current_);
		}
		current_actual_ = false;
	}

	virtual ~QueueIterator() {
	}
};


#endif /* QUEUE_ITERATOR_HPP_ */


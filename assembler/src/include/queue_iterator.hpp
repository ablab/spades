//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef QUEUE_ITERATOR_HPP_
#define QUEUE_ITERATOR_HPP_

#include "verify.hpp"
#include <set>
#include <queue>
#include "standard_base.hpp"

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

/*
 * This map works only if
 * 1. an element that was deleted once never appears again
 * 2. Element is still valid for operations like hash and == even after it was deleted
 */
template<typename Key, typename Comparator>
class fast_erasable_priority_queue {
private:
	std::priority_queue<Key, std::vector<Key>, Comparator> storage_;
	unordered_set<Key> deleted_;

	void update() {
		while(!storage_.empty()) {
			Key top = storage_.top();
			if(deleted_.count(top) > 0) {
				storage_.pop();
			} else {
				return;
			}
		}
	}

public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with VERIFY(false) inside it.
	 */
	fast_erasable_priority_queue(const Comparator& comparator = Comparator()) :
		storage_(comparator) {
	}

	template<typename InputIterator>
	fast_erasable_priority_queue(InputIterator begin, InputIterator end,
			const Comparator& comparator = Comparator()) :
		storage_(begin, end, comparator) {
	}

	void pop() {
		VERIFY(!empty());
		Key top_element = top();
		storage_.pop();
		update();
	}

	const Key& top() const {
		VERIFY(!empty());
		return storage_.top();
	}

	void push(const Key& key) {
		VERIFY(deleted_.count(key) == 0);
		storage_.push(key);
	}

	void erase(const Key& key) {
		deleted_.insert(key);
		update();
	}

	bool empty() const {
		return storage_.empty();
	}

	size_t size() const {
		VERIFY(false);
		return storage_.size();
	}

	template <class InputIterator>
	void insert ( InputIterator first, InputIterator last ) {
		for(; first != last; ++first) {
			storage_.push(*first);
		}
	}
};

template<typename ElementId, typename Comparator = std::less<ElementId> >
class QueueIterator {
private:

	struct ReverseComparator {
		Comparator comparator_;
		ReverseComparator(Comparator comparator) : comparator_(comparator) {
		}
		bool operator()(ElementId a, ElementId b) const {
			return comparator_(b, a);
		}
	};

	bool current_actual_;
	bool current_deleted_;
	ElementId current_;
	fast_erasable_priority_queue<ElementId, ReverseComparator> queue_;
protected:

	template<typename InputIterator>
	void insert(InputIterator begin, InputIterator end) {
		queue_.insert(begin, end);
	}

	QueueIterator(const Comparator& comparator = Comparator()) :
		current_actual_(false), current_deleted_(false), queue_(ReverseComparator(comparator)) {
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
		if(!current_actual_) {
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


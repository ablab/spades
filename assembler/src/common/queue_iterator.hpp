#ifndef QUEUE_ITERATOR_HPP_
#define QUEUE_ITERATOR_HPP_

#include <set>

template<typename Key, typename Comparator = std::less<Key> >
class erasable_priority_queue {
private:
	std::set<Key, Comparator> storage_;
public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with assert(false) inside it.
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
		assert(!storage_.empty());
		storage_.erase(storage_.begin());
	}

	const Key& top() const {
		return *(storage_.begin());
	}

	void push(const Key& key) {
		storage_.insert(key);
	}

	bool erase(const Key& key) {
		return storage_.erase(key) > 0;
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
	bool ready;
	erasable_priority_queue<ElementId, Comparator> queue_;
protected:

	template<typename InputIterator>
	void insert(InputIterator begin, InputIterator end) {
		queue_.insert(begin, end);
	}

	QueueIterator(const Comparator& comparator = Comparator()) :
		ready(true), queue_(comparator) {
	}

	template<typename iterator>
	QueueIterator(iterator begin, iterator end,
			const Comparator& comparator = Comparator()) :
		ready(true), queue_(comparator) {
		fillQueue(begin, end);
	}

	void erase(const ElementId& toRemove) {
		if (ready && toRemove == queue_.top()) {
			ready = false;
		}
		queue_.erase(toRemove);
	}

	void push(const ElementId& toAdd) {
		queue_.push(toAdd);
	}

public:
	bool isEnd() const {
		return queue_.empty();
	}

	ElementId operator*() const {
		assert(!queue_.empty());
		assert(ready);
		return queue_.top();
	}

	void operator++() {
		if (ready) {
			queue_.pop();
		}
		else {
			ready = true;
		}
	}

	virtual ~QueueIterator() {
	}
};


#endif /* QUEUE_ITERATOR_HPP_ */

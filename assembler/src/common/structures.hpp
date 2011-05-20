#ifndef STRUCTURES_HPP_
#define STRUCTURES_HPP_

#include <set>

template<typename Key, typename Comparator = std::less<Key> >
class erasable_priprity_queue {
private:
	std::set<Key, Comparator> storage_;
public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with assert(false) inside it.
	 */
	erasable_priprity_queue(const Comparator& comparator = Comparator())  __attribute__ ((deprecated)) :
		storage_(comparator) {
	}

	template<typename InputIterator>
	erasable_priprity_queue(InputIterator begin, InputIterator end,
			const Comparator& comparator = Comparator()) :
		storage_(begin, end, comparator) {
	}

	Key pop() {
		Key key = top();
		storage_.erase(storage_.begin());
		return key;
	}
	const Key& top() const {
		return *(storage_.begin());
	}

	void push(const Key& key) {
		storage_.insert(key);
	}

	bool remove(const Key& key) {
		return storage_.erase(key) > 0;
	}

	bool empty() const {
		return storage_.empty();
	}

	size_t size() const {
		return storage_.size();
	}
};

// use std::priority_queue iterator !!!
// http://www.cplusplus.com/reference/stl/priority_queue/


template<typename ElementId, typename Comparator = std::less<ElementId> >
class QueueIterator {
private:
	bool ready;
protected:
	erasable_priprity_queue<ElementId, Comparator> queue_;

	template<typename iterator>
	void AddAll(iterator begin, iterator end) {
		for (iterator it = begin; it != end; ++it) {
			queue_.push(*it);
		}
	}

	QueueIterator(const Comparator& comparator = Comparator())  __attribute__ ((deprecated)) :
		ready(true), queue_(comparator) {
	}

	template<typename iterator>
	QueueIterator(iterator begin, iterator end,
			const Comparator& comparator = Comparator()) :
		ready(true), queue_(comparator) {
		fillQueue(begin, end);
	}

	void remove(ElementId toRemove) {
		if (ready && toRemove == queue_.top()) {
			ready = false;
		}
		queue_.remove(toRemove);
	}

public:
	//== is supported only in case this or other is end iterator
	bool operator==(const QueueIterator& other) {
		if (this->queue_.empty() && other.queue_.empty())
			return true;
		if (this->queue_.empty() || other.queue_.empty())
			return false;
		assert(false);
	}

	bool operator!=(const QueueIterator& other) {
		if (this->queue_.empty() && other.queue_.empty())
			return false;
		if (this->queue_.empty() || other.queue_.empty())
			return true;
		assert(false);
	}

	ElementId operator*() const {
		assert(!queue_.empty());
		assert(ready);
		return queue_.top();
	}

	void operator++() {
		assert(!queue_.empty());
		if (ready) {
			queue_.pop();
		}
		else {
			ready = true;
		}
		//		cout << "remove " << queue_.size() << endl;
	}

	virtual ~QueueIterator() {
	}
};


#endif /* STRUCTURES_HPP_ */

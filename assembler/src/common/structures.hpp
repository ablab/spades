#ifndef STRUCTURES_HPP_
#define STRUCTURES_HPP_

// use std::priority_queue !!!
// http://www.cplusplus.com/reference/stl/priority_queue/

template<typename Key, typename Comparator = std::less<Key> >
class PriorityQueue {
private:
	set<Key, Comparator> storage_;
public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with assert(false) inside it.
	 */
	PriorityQueue(const Comparator& comparator = Comparator())  __attribute__ ((deprecated)) :
		storage_(comparator) {
	}

	template<typename InputIterator>
	PriorityQueue(InputIterator begin, InputIterator end,
			const Comparator& comparator = Comparator()) :
		storage_(begin, end, comparator) {
	}

	Key poll() {
		Key key = *(storage_.begin());
		storage_.erase(storage_.begin());
		return key;
	}
	Key peek() const {
		return *(storage_.begin());
	}

	void offer(const Key key) {
		storage_.insert(key);
	}

	bool remove(const Key key) {
		return storage_.erase(key) > 0;
	}

	bool empty() const {
		return storage_.empty();
	}

	size_t size() const {
		return storage_.size();
	}
};

template<typename ElementId, typename Comparator = std::less<ElementId> >
class QueueIterator {
private:
	bool ready;
protected:
	PriorityQueue<ElementId, Comparator> queue_;

	template<typename iterator>
	void AddAll(iterator begin, iterator end) {
		for (iterator it = begin; it != end; ++it) {
			queue_.offer(*it);
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
		if (ready && toRemove == queue_.peek()) {
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
		return queue_.peek();
	}

	void operator++() {
		assert(!queue_.empty());
		if (ready)
			queue_.poll();
		else
			ready = true;
		//		cout << "remove " << queue_.size() << endl;
	}

	virtual ~QueueIterator() {
	}
};


#endif /* STRUCTURES_HPP_ */

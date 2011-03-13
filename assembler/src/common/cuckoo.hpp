/*
 * cuckoo.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 */

#ifndef CUCKOO_HPP_
#define CUCKOO_HPP_

using namespace std;

template <class Key, class Value, class Hash1, class Hash2, class Hash3, class Pred>
class cuckoo {
private:
	struct Data {
		pair<Key, Value> p;
		bool exists;
		Data() : exists(false) {}
	};
	Data *data_;
	size_t len_;
	size_t size_;
public:
	class iterator {
	private:
		size_t pos;
		cuckoo *hash;
		iterator(size_t p, cuckoo *h) : pos(p), hash(h) {}
		friend class cuckoo;
	public:
		iterator() : pos(0), hash(NULL) {
			// nothing
		}
		void operator=(const iterator &it) {
			pos = it.pos;
			hash = it.hash;
		}
		iterator& operator++() {
			assert(hash != NULL);
			assert(pos != hash->len_);
			while (pos < len_ && !data_[pos].exists) {
				pos++;
			}
			return *this;
		}
		iterator operator++(int) {
			iterator res = *this;
			this->operator++();
			return res;
		}
		pair<Key, Value>* operator->() {
			assert(pos != hash->len_);
			return &hash->data_[pos].p;
		}
		bool operator==(const iterator &it) {
			return pos == it.pos && hash == it.hash;
		}
		bool operator!=(const iterator &it) {
			return !(*this == it);
		}
	};
private:

	bool is_here(const Key &k, size_t pos) {
		return data_[pos].exists && Pred()(data_[pos].p.first, k);
	}

	// does it needs optimization (for not creating Hash & Pred objects)?
	size_t hash(const Key &k, size_t hash_num) {
		if (hash_num == 1) {
			return Hash1()(k) % size_;
		}
		if (hash_num == 2) {
			return Hash2()(k) % size_;
		}
		if (hash_num == 3) {
			return Hash3()(k) % size_;
		}
		assert(false); // should never happens
		return -1;
	}

	/*pair<const Key,Value>& pair<const Key, Value>::operator=(const pair<const Key, Value> &p) {

	}*/

	// this one should increase len_ and reallocate data_ is necessary!
	iterator add_new(pair<Key, Value> p) {
		size_t blocked_pos = -1;
		while (true) {
			for (size_t i = 1; i <= 3; ++i) {
				size_t pos = hash(p.first, i);
				if (pos == blocked_pos) continue;
				if (!data_[pos].exists) {
					data_[pos].p = p;
					data_[pos].exists = true;
					return iterator(pos, this);
				}
			}
			// no free space, put it in random position and repeat with another element
			size_t pos;
			do {
				pos = rand() % 3 + 1;
			} while (pos == blocked_pos);
			assert(data_[pos].exists);
			pair<const Key, Value> p_next = data_[pos].p;
			data_[pos].p = pair<const Key, Value>(p.first, p.second);
			p = p_next;
			blocked_pos = pos;
		}
		assert(false); // should never happens
		return end();
	}

public:
	cuckoo(size_t len) {
		data_ = new Data[len]; // malloc?
	}

	~cuckoo() {
		delete[] data_;
	}

	iterator begin() {
		return iterator(0);
	}

	iterator end() {
		return iterator(len_, this);
	}

	iterator find(const Key &k) {
		for (size_t i = 1; i <= 3; ++i) {
			size_t pos = hash(k, i);
			if (is_here(k, pos)) {
				return iterator(pos, this);
			}
		}
		return end();
	}

	// TODO: increase data len and rehash if necessary!
	pair<iterator, bool> insert(const pair<const Key, Value> &k) {
		iterator res = find(k.first);
		if (res != end()) {
			res->second = k.second;
			return make_pair(res, false);
		}
		assert(res == end());
		res = add_new(k);
		return make_pair(res, false);
	}

	size_t size() const {
		return size_;
	}

};



#endif /* CUCKOO_HPP_ */

/*
 * cuckoo.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 *  Last modify: 14.02.2011 00:35
 *      Author: Mariya Fomkina
 */

#ifndef CUCKOO_HPP_
#define CUCKOO_HPP_

using namespace std;

// Classes Key and Value are some data types, that should be stored in hash
// (e.g. string and int). 
// Class Hash is a function that gets data of type Key and some number 
// (which is the number of appropriate hash function) and returns size_t 
// (e.g. Hash(int, size_t)).
// Class Pred is predicator that compares two Key values.
// Value size_t d is the number of hash functions (thus arrays also) 
// that will be used in the program (can be >= 2).
// Value size_t init_length is the initial length of the whole structure.
// Value size_t max_loop determines the maximum number of kick cycles during 
// insertion before rehash.
// Example of use: 
// cuckoo<int, int, Hasher, std::equal_to<int>, 4, 10, 100> Cuckoo; 
template <class Key, class Value, class Hash, class Pred, size_t d = 3, 
	  size_t init_length = 15, size_t max_loop = 100>
class cuckoo {
private:
  struct Data {
    pair<Key, Value> p;
    bool exists;
    Data() : exists(false) {}
  };
  // The array of vectors, each of which is hash array 
  vector<Data> *data_[d];
  // The total length of all the hash arrays
  size_t len_;
  // The length of every hash array
  size_t len_part_;
  // The actual number of elements in cuckoo hash
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
      ++pos;
      while (pos < hash->len_ && !hash->data_from(pos).exists) {
	++pos;
      }
      return *this;
    }
    iterator operator++(int) {
      iterator res = *this;
      this->operator++();
      return res;
    }
    pair<Key, Value>& operator*() {
      return (*hash).data_from(pos).p;
    }
    bool operator==(const iterator &it) {
      return pos == it.pos && hash == it.hash;
    }
    bool operator!=(const iterator &it) {
      return !(*this == it);
    }
  };

private:

  Data& data_from(size_t pos) {
    return (*(data_[pos / len_part_]))[pos % len_part_];
  }
  
  bool is_here(const Key &k, size_t pos) {
    return data_from(pos).exists && Pred()(data_from(pos).p.first, k);
  }
  
  size_t hash(const Key &k, size_t hash_num) {
    return Hash()(k, hash_num) % len_part_;
  }
  
  /*pair<const Key,Value>& pair<const Key, Value>::operator=(const pair<const Key, Value> &p) {
    }*/

  // Works incorrectly, contain debug output.
  // I will try to fix bug as soon as possible.
  void rehash() {
    cerr << "Started rehashing with " << size () << 
      " and " << size_ << endl;
    len_part_ = len_part_ * 6 / 5 + 1;
    len_ = len_part_ * d;
    for (size_t i = 0; i < d; ++i) {
      (*(data_[i])).resize(len_part_);
    } 
    iterator it = begin();
    if (!data_from(it.pos).exists) ++it;
    while (it != end()) {
      size_t i = it.pos / len_part_;
      size_t j = it.pos % len_part_;
      if (j != hash((*it).first, i)) {
	pair<Key, Value> t = *it;
	remove(it);
	add_new(t);
      } else { 
	++it;
      }
    }
    cerr << "    and finished with " << size() << endl;
  }

  iterator add_new(pair<Key, Value> p) {
    for (size_t i = 0; i < max_loop; ++i) {
      for (size_t j = 0; j < d; ++j) {
	size_t pos = hash(p.first, j);
	swap(p, data_from(j * len_part_ + pos).p);
	bool exists = data_from(j * len_part_ + pos).exists;
	data_from(j * len_part_ + pos).exists = true;
	if (!exists) {
	  return iterator(j * len_part_ + pos, this);
	} 
      }
    }
    rehash();
    add_new(p);
    return end(); 
  }
  
public:
  cuckoo() {
    len_part_ = init_length / d + 1;
    len_ = len_part_ * d;
    for (size_t i = 0; i < d; ++i) {
      data_[i] = new vector<Data>(len_part_);
    }
    size_ = 0;
  }
  
  ~cuckoo() {
    for (size_t i = 0; i < d; ++i) {
      delete data_[i];
    }
  }

  iterator begin() {
    return iterator(0, this);
  }
  
  iterator end() {
    return iterator(len_, this);
  }
  
  iterator remove(iterator & it) {
    data_from(it.pos).exists = false;
    return ++it;
  }

  iterator find(const Key &k) {
    for (size_t i = 0; i < d; ++i) {
      size_t pos = hash(k, i);
      if (is_here(k, i * len_part_ + pos)) {
	return iterator(i * len_part_ + pos, this);
      }
    }
    return end();
  }
  
  // Returns iterator to the value and true if new value was inserted
  // and false otherwise.
  pair<iterator, bool> insert(const pair<const Key, Value> &k) {
    iterator res = find(k.first);
    if (res != end()) {
      (*res).second = k.second;
      return make_pair(res, false);
      } 
    assert(res == end());
    res = add_new(k);
    ++size_;
    return make_pair(res, true);
  }
  
  // temporary variant of function - in debug it is
  // possible to see that smth is wrong with size
  size_t size() {
    size_t k = 0;
    for (size_t i = 0; i < len_; ++i) 
      if (data_from(i).exists) ++k; 
    return k;
    //return size_;
  }

  size_t length() const {
    return len_;
  }
};

#endif /* CUCKOO_HPP_ */

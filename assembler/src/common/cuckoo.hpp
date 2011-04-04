/*
 * cuckoo.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 *  Last modify: 28.03.2011 21:09
 *      Author: Mariya Fomkina
 */

#include <cstring>

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
// When you know the approximate number of records to be used, 
// it is a good idea to take this value in 1.05-1.1 times more and
// small value of step.
// Value size_t max_loop determines the maximum number of kick cycles during 
// insertion before rehash.
// Values size_t step_nom and step_denom determine the ratio of 
// increasing the size of hash during rehash (step_nom/step_denom). 
// The less it is the less memory will be used but the more time is needed. 
// Example of use: 
// cuckoo<int, int, Hasher, std::equal_to<int>, 4, 10, 50, 6, 5> Cuckoo; 
template <class Key, class Value, class Hash, class Pred, size_t d = 4, 
	  size_t init_length = 100, size_t max_loop = 100, 
	  size_t step_nom = 3, size_t step_denom = 2>
class cuckoo {
private:
  // The array of vectors, each of which is hash array 
  typedef pair<Key, Value> Data;
  Data** data_;
  // The array of flags indicating existence of the element in hash
  char* exists_;
  // The total length of all the hash arrays
  size_t len_;
  // The length of every hash array
  size_t len_part_;
  // The actual number of elements in cuckoo hash
  size_t size_;
  // The flag that anounces that rehash was made recently
  bool is_rehashed_;
public:
  friend class iterator;

  class iterator {
  private:
    size_t pos;
    cuckoo* hash;
    iterator(size_t p, cuckoo* h): pos(p), hash(h) {}
    friend class cuckoo;
  public:
    iterator() : pos(0), hash(NULL) {}

    void operator=(const iterator &it) {
      pos = it.pos;
      hash = it.hash;
    }

    iterator& operator++() {
      assert(hash != NULL);
      assert(pos != hash->len_);
      ++pos;
      while (!hash->get_exists(pos) && pos < hash->len_) {
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
      return (*hash).data_from(pos);
    }

    bool operator==(const iterator &it) {
      return pos == it.pos && hash == it.hash;
    }

    bool operator!=(const iterator &it) {
      return !(*this == it);
    }
  };

private:

  bool get_exists(size_t pos) const {
    return (bool)(exists_[pos / 8] & (1 << (7 - (pos % 8))));  
  }

  void set_exists(size_t pos) {
    exists_[pos / 8] = exists_[pos / 8] | (1 << (7 - (pos % 8)));
  }

  void unset_exists(size_t pos) {
    exists_[pos / 8] = exists_[pos / 8] & (~(1 << (7 - (pos % 8))));
  }

  void init() {
    len_part_ = init_length / d + 1;
    len_part_ = ((len_part_ + 7) / 8) * 8;
    len_ = len_part_ * d;
    data_ = new Data*[d];
    for (size_t i = 0; i < d; ++i) {
      data_[i] = new Data[len_part_];
    }
    exists_ = new char[len_ / 8]; 
    for (size_t i = 0; i < len_ / 8; ++i) exists_[i] = 0;
    size_ = 0;
    is_rehashed_ = false;
  }

  void clear_all() {
    for (size_t i = 0; i < d; ++i) {
      delete [] data_[i];
    }
    delete [] data_;
    delete [] exists_;
  }

  inline Data& data_from(size_t pos) const {
    return data_[pos / len_part_][pos % len_part_];
  }
  
  inline bool is_here(const Key &k, size_t pos) const {
    return get_exists(pos) && Pred()(data_from(pos).first, k); 
  }

  inline size_t hash(const Key &k, size_t hash_num) {
    return Hash()(k, hash_num) % len_part_;
  }
  
  void update_exists(size_t len_temp_) { 
    char* t = new char[len_ / 8];
    char* s = t;
    char* f = s + len_ / 8;
    for (; s < f; ++s) {
      *s = 0;
    } 
    for (size_t i = 0; i < d; ++i) {
      memcpy(t + (i * len_part_ / 8), exists_ + (i * len_temp_ / 8), (len_temp_ / 8));
    }
    swap(t, exists_);
    delete [] t;
  }

  void update_data(size_t len_temp_) {
    for (size_t i = 0; i < d; ++i) {
      Data* t = new Data[len_part_];
      memcpy(t, data_[i], len_temp_*sizeof(Data));
      swap(t, data_[i]);
      delete [] t;
    } 
  }

  void rehash() {
    size_t len_temp_ = len_part_;
    len_part_ = len_part_ * step_nom / step_denom;
    len_part_ = ((len_part_ + 7) / 8) * 8;
    len_ = len_part_ * d;
    
    update_exists(len_temp_);
    update_data(len_temp_);

    iterator it = begin();
    if (!get_exists(it.pos)) ++it;
    while (it != end()) {
      size_t i = it.pos / len_part_;
      size_t j = it.pos % len_part_;
      if (j != hash((*it).first, i)) {
	Data t = *it;
	remove(it);
	add_new(t);
	if (is_rehashed_) {
	  it = begin();
	  if (!get_exists(it.pos)) ++it; 
	}
      } else { 
	++it;
      }
    }
    is_rehashed_ = true;
  }

  iterator add_new(pair<Key, Value> p) {
    for (size_t i = 0; i < max_loop; ++i) {
      for (size_t j = 0; j < d; ++j) {
	size_t pos = hash(p.first, j);
	swap(p, data_from(j * len_part_ + pos));
	bool exists = get_exists(j * len_part_ + pos); 
	set_exists(j * len_part_ + pos);
	if (!exists) {
	  is_rehashed_ = false;
	  ++size_;
	  return iterator(j * len_part_ + pos, this);
	} 
      }
    }
    rehash();
    return add_new(p);
  }

  iterator remove(iterator& it) {
    unset_exists(it.pos);
    --size_;    
    return ++it;
  }

public:
  cuckoo() {
    init();
  }
  
  ~cuckoo() {
    clear_all();
  }

  cuckoo<Key, Value, Hash, Pred, d, init_length, max_loop, step_nom, step_denom>& operator=
  (cuckoo<Key, Value, Hash, Pred, d, init_length, max_loop, step_nom, step_denom>& Cuckoo) {
    clear_all();
    init();
    iterator it = Cuckoo.begin();
    if (!(Cuckoo.get_exists(it.pos))) ++it;
    iterator final = Cuckoo.end();
    while (it != final) {
      insert(*it);
      ++it;
      }
    return *this;
  }

  cuckoo(cuckoo<Key, Value, Hash, Pred, d, init_length, max_loop, step_nom, step_denom>& Cuckoo) {
    init();
    *this = Cuckoo;
  }

  inline iterator begin() {
    return iterator(0, this);
  }
  
  inline iterator end() {
    return iterator(len_, this);
  }

  Value& operator[](const Key& k) {
    iterator it = find(k);
    if (it == end()) {
      it = insert(make_pair(k, Value())).first;
    }
    return (*it).second;
  }

  void erase(iterator& it) {
    remove(it);
  }

  void erase(iterator& first, iterator& last) {
    while (first != last) {
      first = remove(first);
    }
  }

  size_t erase(const Key& k) {
    iterator it = find(k);
    if (it != end()) {
      remove(it);
      return 1;
    }
    return 0;
  }

  iterator find(const Key& k) {
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
    return make_pair(res, true);
  }

  void clear() {
    char* t = new char[len_ / 8];
    for (size_t i = 0; i < len_ / 8; ++i) t[i] = 0;
    swap(t, exists_);
    size_ = 0;
  }

  inline bool empty() {
    return (size_ == 0);
  }
 
  inline size_t size() const {
    return size_;
  }

  inline size_t length() const {
    return len_;
  }
};

#endif /* CUCKOO_HPP_ */

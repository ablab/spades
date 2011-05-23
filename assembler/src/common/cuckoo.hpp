/*
 * cuckoo.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 *  Last modify: 15.05.2011
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
// Example of use: 
// cuckoo<int, int, Hasher, std::equal_to<int> > Cuckoo; 
template <class Key, class Value, class Hash, class Pred> 
class cuckoo {
private:
  // Next 4 parameters described near constructor
  size_t d_;
  size_t init_length_;
  size_t max_loop_;
  double step_;
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
  friend class const_iterator;

  class const_iterator {
  private:
    size_t pos;
    const cuckoo* hash;
    const_iterator(const size_t p, const cuckoo* h) : pos(p), hash(h) {}
    friend class cuckoo;

  public:
    const_iterator() : pos(0), hash(NULL) {}
    
    void operator=(const const_iterator &it) {
      pos = it.pos;
      hash = it.hash;
    }

    const_iterator& operator++() {
      assert(hash != NULL);
      assert(pos != hash->len_);
      ++pos;
      while ((pos < hash->len_) && (!(hash->get_exists(pos)))) {
        ++pos;
      }
      return *this;
    }
    
    const_iterator operator++(int) {
      const_iterator res = *this;
      this->operator++();
      return res;
    } 

    const pair<Key, Value>& operator*() const {
      return (*hash).data_from(pos);
    }

    const pair<Key, Value>* operator->() const {
      return &((*hash).data_from(pos));
    }

    bool operator==(const const_iterator &it) {
      return pos == it.pos && hash == it.hash;
    }

    bool operator!=(const const_iterator &it) {
      return !(*this == it);
    }
  };

  class iterator {
  private:
    size_t pos;
    cuckoo* hash;
    iterator(size_t p, cuckoo* h) : pos(p), hash(h) {}
    friend class cuckoo;

  public:
    iterator() : pos(0), hash(NULL) {}

    operator const_iterator() {
      return const_iterator(pos, hash);
    }

    void operator=(const iterator &it) {
      pos = it.pos;
      hash = it.hash;
    }

    iterator& operator++() {
      assert(hash != NULL);
      assert(pos != hash->len_);
      ++pos;
      while ((pos < hash->len_) && !(hash->get_exists(pos))) {
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

    pair<Key, Value>* operator->() {
      return &((*hash).data_from(pos));
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
    len_part_ = init_length_ / d_ + 1;
    len_part_ = ((len_part_ + 7) / 8) * 8;
    len_ = len_part_ * d_;
    data_ = new Data*[d_];
    for (size_t i = 0; i < d_; ++i) {
      data_[i] = new Data[len_part_];
    }
    exists_ = new char[len_ / 8]; 
    for (size_t i = 0; i < len_ / 8; ++i) exists_[i] = 0;
    size_ = 0;
    is_rehashed_ = false;
  }

  void clear_all() {
    for (size_t i = 0; i < d_; ++i) {
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

  inline size_t hash(const Key &k, size_t hash_num) const {
    return Hash()(k, hash_num) % len_part_;
  }
  
  void update_exists(size_t len_temp_) { 
    char* t = new char[len_ / 8];
    char* s = t;
    char* f = s + len_ / 8;
    for (; s < f; ++s) {
      *s = 0;
    } 
    for (size_t i = 0; i < d_; ++i) {
      memcpy(t + (i * len_part_ / 8), exists_ + (i * len_temp_ / 8), (len_temp_ / 8));
    }
    std::swap(t, exists_);
    delete [] t;
  }

  void update_data(size_t len_temp_) {
    for (size_t i = 0; i < d_; ++i) {
      Data* t = new Data[len_part_];
      memcpy(t, data_[i], len_temp_*sizeof(Data));
      std::swap(t, data_[i]);
      delete [] t;
    } 
  }

  void rehash() {
    size_t len_temp_ = len_part_;

    len_part_ = (size_t)(len_part_ * step_);
    len_part_ = ((len_part_ + 7) / 8) * 8;
    len_ = len_part_ * d_;
    
    update_exists(len_temp_);
    update_data(len_temp_);

    size_t n = 0;
    iterator it = begin();
    while (it != end()) {
      size_t i = it.pos / len_part_;
      size_t j = it.pos % len_part_;
      ++n;
      if (j != hash((*it).first, i)) {
        Data t = *it;
        remove(it);
        add_new(t);
        if (is_rehashed_) {
          it = begin();
          is_rehashed_ = false;
        }
      } else { 
        ++it;
      }
    }
    is_rehashed_ = true;
  }

  iterator add_new(pair<Key, Value> p) {
    for (size_t i = 0; i < max_loop_; ++i) {
      for (size_t j = 0; j < d_; ++j) {
        size_t pos = hash(p.first, j);
        std::swap(p, data_from(j * len_part_ + pos));
        bool exists = get_exists(j * len_part_ + pos); 
        set_exists(j * len_part_ + pos);
        if (!exists) {
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
  // @parameter d the number of hash functions (thus arrays also) 
  // that will be used in the program (can be >= 2).
  // @parameter init_length the initial length of the whole structure.
  // When you know the approximate number of records to be used, 
  // it is a good idea to take this value in 1.05-1.1 times more and
  // small value of step.
  // @parameter max_loop determines the maximum number of kick cycles during 
  // insertion before rehash.
  // @parameter step determines the ratio of increasing the size of hash
  // during rehash.   
  // The less it is the less memory will be used but the more time is needed. 
  cuckoo(size_t d = 4, size_t init_length = 100, 
         size_t max_loop = 100, double step = 1.2)
    : d_(d), init_length_(init_length), 
      max_loop_(max_loop), step_(step) {
    init();
  }
  
  ~cuckoo() {
    clear_all();
  }

  cuckoo<Key, Value, Hash, Pred>& operator=
  (const cuckoo<Key, Value, Hash, Pred>& Cuckoo) {
    clear_all();
    d_ = Cuckoo.d_;
    init_length_ = Cuckoo.init_length_;
    max_loop_ = Cuckoo.max_loop_;
    step_ = Cuckoo.step_;
    init();
    const_iterator it = Cuckoo.begin();
    const_iterator final = Cuckoo.end();
    while (it != final) {
      insert(*it);
      ++it;
    }
    return *this;
  }

  //TODO: optimize!
  cuckoo(const cuckoo<Key, Value, Hash, Pred>& Cuckoo) {
    init();
    *this = Cuckoo;
  }

  template <class InputIterator>
  cuckoo(InputIterator first, InputIterator last) {
    d_ = 4;
    init_length_ = 100;
    max_loop_ = 100;
    step_ = 1.2;
    init();
    for (iterator it = first; it != last; ++it) {
      add_new(*it);
    }
  }

  // For test only!!!
  void set_up(size_t d = 4, size_t init_length = 100, 
              size_t max_loop = 100, double step = 1.2) {
    clear_all();
    d_ = d;
    init_length_ = init_length;
    max_loop_ = max_loop;
    step_ = step;
    init();
  }

  void swap(cuckoo<Key, Value, Hash, Pred> Cuckoo) {
    std::swap(d_, Cuckoo.d_);
    std::swap(init_length_, Cuckoo.init_length_);
    std::swap(max_loop_, Cuckoo.max_loop_);
    std::swap(step_, Cuckoo.step_);
    std::swap(data_, Cuckoo.data_);
    std::swap(exists_, Cuckoo.exists_);
    std::swap(len_, Cuckoo.len_);
    std::swap(len_part_, Cuckoo.len_part_);
    std::swap(size_, Cuckoo.size_);
  }

  inline iterator begin() {
    iterator it = iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }
  
  inline const_iterator begin() const {
    const_iterator it = const_iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }
  
  inline iterator end() {
    return iterator(len_, this);
  }

  inline const_iterator end() const {
    return const_iterator(len_, this);
  }

  Value& operator[](const Key& k) {
    iterator it = find(k);
    if (it == end()) {
      it = insert(make_pair(k, Value())).first;
    }
    return (*it).second;
  }

  void erase(iterator it) {
    remove(it);
  }

  void erase(iterator first, iterator last) {
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
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      if (is_here(k, i * len_part_ + pos)) {
        return iterator(i * len_part_ + pos, this);
      }
    }
    return end();
  }
  
  const_iterator find(const Key& k) const {
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      if (is_here(k, i * len_part_ + pos)) {
        return const_iterator(i * len_part_ + pos, this);
      }
    }
    return end();
  }

  size_t count(const Key& k) const {
    if (find(k) != end()) {
      return 1;
    } else {
      return 0;
    }
  }
  
  pair<iterator, iterator> equal_range(const Key& k) {
    if (find(k) != end()) {
      return std::make_pair<iterator, iterator>(find(k), find(k));
    } else {
      return std::make_pair<iterator, iterator>(end(), end());
    }
  }

  pair<const_iterator, const_iterator> equal_range(const Key& k) const {
    if (find(k) != end()) {
      return std::make_pair<const_iterator, const_iterator>
        (find(k), find(k));
    } else {
      return std::make_pair<const_iterator, const_iterator>
        (end(), end());
    }
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

  template <class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    for (iterator it = first; it != last; ++it) {
      insert(*it);
    }
  } 

  void clear() {
    //char* t = new char[len_ / 8];
    //for (size_t i = 0; i < len_ / 8; ++i) t[i] = 0;
    //swap(t, exists_);
    //size_ = 0;
    clear_all();
    init();
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

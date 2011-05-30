/**
 * @file    cuckoo.hpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Cuckoo hashing implementation with the interface close to
 * standard std::map and std::unordered_map interface.
 *
 * See http://en.wikipedia.org/wiki/Cuckoo_hashing
 */

#include <cstring>

#ifndef _CUCKOO_HPP_
#define _CUCKOO_HPP_

const static size_t D = 3;
const static size_t INIT_LENGTH = 100;
const static size_t MAX_LOOP = 100;
const static double STEP = 1.2;

/**
 * @param Key key type in hash
 * @param Value value type in hash
 * @param Hash function that gets data of type Key and some number 
 * (which is the number of appropriate hash function) and returns size_t 
 * (e.g. Hash(int, size_t))
 * @param Equal predicator that compares two Key values
 * @example cuckoo<int, int, Hasher, std::key_equal_to<int> > Cuckoo;
 */
template <class Key, class Value, class Hash, class Equal> 
class cuckoo {
private:
  typedef pair<Key, Value> Data;

  size_t d_;
  size_t init_length_;
  size_t max_loop_;
  double step_;
  Hash hasher_;
  Equal key_equal_;

  /** 
   * @variable The array of vectors, each of which is hash array
   */ 
  Data** data_;

  /**
   * @variable The array of flags indicating existence of the element in hash
   */  
  char* exists_;

  /**
   * @variable The total length of all the hash arrays
   */
  size_t len_;

  /**
   * @variable The length of every hash array
   */
  size_t len_part_;

  /**
   * @variable The actual number of elements in cuckoo hash
   */  
  size_t size_;

  /**
   * @variable The flag that anounces that rehash was made recently
   */
  bool is_rehashed_;

public:
  friend class iterator;
  friend class const_iterator;

  /**
   * Used with const cuckoo objects.
   */
  class const_iterator {
  private:
    size_t pos;
    const cuckoo* hash;
    const_iterator(const size_t p, const cuckoo* h) : pos(p), hash(h) {}

  public:
    friend class cuckoo;

    /**
     * Default constructor.
     */
    const_iterator() : pos(0), hash(NULL) {}

    void operator=(const const_iterator& it) {
      pos = it.pos; 
      hash = it.hash;
    }

    const_iterator(const const_iterator& it) {
      *this = it;
    }

    const_iterator& operator++() {
      assert(hash != NULL);
      if (pos != hash->len_) ++pos;
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

    const Data& operator*() const {
      return (*hash).data_from(pos);
    }

    const Data* operator->() const {
      return &((*hash).data_from(pos));
    }

    bool operator==(const const_iterator &it) {
      return (pos == it.pos) && (hash == it.hash);
    }

    bool operator!=(const const_iterator &it) {
      return !(*this == it);
    }
  };

  /**
   * Used with non-const cuckoo objects.
   */
  class iterator {
  private:
    size_t pos;
    cuckoo* hash;
    iterator(size_t p, cuckoo* h) : pos(p), hash(h) {}

  public:
    friend class cuckoo;

    /**
     * Default constructor.
     */
    iterator() : pos(0), hash(NULL) {}

    /**
     * Convertion to const_iterator.
     */
    operator const_iterator() {
      return const_iterator(pos, hash);
    }

    void operator=(const iterator &it) {
      pos = it.pos;
      hash = it.hash;
    }

    iterator(const iterator &it) {
      *this = it;
    }

    iterator& operator++() {
      assert(hash != NULL);
      if (pos != hash->len_) ++pos;
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

    Data& operator*() {
      return (*hash).data_from(pos);
    }

    Data* operator->() {
      return &((*hash).data_from(pos));
    }

    bool operator==(const iterator &it) {
      return (pos == it.pos) && (hash == it.hash);
    }

    bool operator!=(const iterator &it) {
      return !(*this == it);
    }
  };

private:

  /**
   * Check whether there is a Data element at pos.
   *
   * @param pos Position in hash arrays
   */
  bool get_exists(size_t pos) const {
    return (bool)(exists_[pos >> 3] & (1 << (pos - ((pos >> 3) << 3))));
  }

  /**
   * Set flag of existence of Data element at pos.
   *
   * @param pos Position in hash arrays
   */
  void set_exists(size_t pos) {
    exists_[pos >> 3] = exists_[pos >> 3] | (1 << (pos - ((pos >> 3) << 3)));
  }

  /**
   * Unset flag of existence Data element at pos.
   *
   * @param pos Position in hash arrays
   */
  void unset_exists(size_t pos) {
    exists_[pos >> 3] = exists_[pos >> 3] & (~(1 << (pos - ((pos >> 3) << 3))));
  }

  /**
   * Initialize all the variables of cuckoo.
   */
  void init() {
    len_part_ = init_length_ / d_ + 1;
    len_part_ = ((len_part_ + 7) >> 3) << 3;
    len_ = len_part_ * d_;
    data_ = (Data**)(malloc(d_ * sizeof(Data*)));
    for (size_t i = 0; i < d_; ++i) {
      data_[i] = (Data*)(malloc(len_part_ * sizeof(Data)));
    }
    exists_ = (char*)(malloc((len_ >> 3) * sizeof(char)));
    for (size_t i = 0; i < len_ >> 3; ++i) exists_[i] = 0;
    size_ = 0;
    is_rehashed_ = false;
  }

  /**
   * Copies information from Cuckoo to clean object.
   */
  void copy(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    d_ = Cuckoo.d_;
    init_length_ = Cuckoo.init_length_;
    max_loop_ = Cuckoo.max_loop_;
    step_ = Cuckoo.step_;
    hasher_ = Cuckoo.hasher_;
    key_equal_ = Cuckoo.key_equal_;
    init();
    const_iterator it = Cuckoo.begin();
    size_t final_len = Cuckoo.length();
    while (it.pos != final_len) {
      insert(*it);
      ++it;
    }
  }

  /**
   * Clear all the data from cuckoo.
   */
  void clear_all() {
    for (size_t i = 0; i < d_; ++i) {
      free(data_[i]);
    }
    free(data_);
    free(exists_);
  }

  /**
   * Get Data element from pos.
   *
   * @param pos Position in hash arrays
   */
  inline Data& data_from(size_t pos) const {
    return data_[pos / len_part_][pos % len_part_];
  }

  /**
   * Check whether element at pos has key k.
   *
   * @param k Key value
   * @param pos Position in hash arrays
   */
  inline bool is_here(const Key& k, size_t pos) const {
    return get_exists(pos) && key_equal_(data_from(pos).first, k);
  }

  /**
   * Return hash function result for key k.
   *
   * @param k Key value
   * @param hash_num The number of hash function from Hash family
   */
  inline size_t hash(const Key& k, size_t hash_num) const {
    return hasher_(k, hash_num) % len_part_;
  }
  
  /**
   * Increase size of exists_ up to len_temp_.
   *
   * @param len_temp_ New size of exists_
   */
  void update_exists(size_t len_temp_) { 
    char* t = (char*)(malloc((len_ >> 3) * sizeof(char)));
    char* s = t;
    char* f = s + (len_ >> 3);
    for (; s < f; ++s) {
      *s = 0;
    } 
    for (size_t i = 0; i < d_; ++i) {
      memcpy(t + (i * (len_part_ >> 3)), exists_ + (i * (len_temp_ >> 3)), (len_temp_ >> 3));
    }
    std::swap(t, exists_);
    free(t);
  }

  /**
   * Increase size of data_ up to len_temp_.
   *
   * @param len_temp_ New size of data_
   */
  void update_data(size_t len_temp_) {
    for (size_t i = 0; i < d_; ++i) {
      Data* t = (Data*)(malloc(len_part_ * sizeof(Data)));
      memcpy(t, data_[i], len_temp_*sizeof(Data));
      std::swap(t, data_[i]);
      free(t);
    } 
  }

  /**
   * Rehash all the cuckoo (i.e. change size and replace Data elements).
   */
  void rehash() {
    size_t len_temp_ = len_part_;

    len_part_ = (size_t)(len_part_ * step_);
    len_part_ = ((len_part_ + 7) >> 3) << 3;
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

  /**
   * Add new Data element.
   *
   * @param p New element
   */
  size_t add_new(Data p) {
    for (size_t i = 0; i < max_loop_; ++i) {
      for (size_t j = 0; j < d_; ++j) {
        size_t pos = hash(p.first, j);
        std::swap(p, data_[j][pos]);
        bool exists = get_exists(j * len_part_ + pos); 
        set_exists(j * len_part_ + pos);
        if (!exists) {
          ++size_;
          return 0;
        } 
      }
    }
    rehash();
    return add_new(p);
  }

  /**
   * Remove element from cuckoo.
   *
   * @param it Iterator to element to be removed
   */
  iterator remove(iterator& it) {
    unset_exists(it.pos);
    --size_;    
    return ++it;
  }

public:
  /** 
   * Default constructor.
   *
   * @param d The number of hash functions (thus arrays also) 
   * that will be used in the program (can be >= 2).
   * @param init_length The initial length of the whole structure.
   * When you know the approximate number of records to be used, 
   * it is a good idea to take this value in 1.05-1.1 times more and
   * small value of step.
   * @param max_loop The maximum number of kick cycles during 
   * insertion before rehash.
   * @param step The ratio of increasing the size of hash during rehash. 
   * @param hasher The hash function object (template parameter by default).
   * @param equal The equal predicator object (template parameter by default).  
   * The less it is the less memory will be used but the more time is needed. 
   */  
  explicit cuckoo(size_t d = D, size_t init_length = INIT_LENGTH, 
                  size_t max_loop = MAX_LOOP, double step = STEP, 
                  const Hash& hasher = Hash(), 
                  const Equal& equal = Equal())
    : d_(d), init_length_(init_length), 
      max_loop_(max_loop), step_(step), 
      hasher_(hasher), key_equal_(equal) {
    init();
  }
  
  /**
   * Destructor.
   */
  ~cuckoo() {
    clear_all();
  }

  /**
   * Copies information from cuckoo of the same type.
   *
   * @param Cuckoo The source of data
   */
  cuckoo<Key, Value, Hash, Equal>& operator=
  (const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    clear_all();
    copy(Cuckoo);
    return *this;
  }

  /**
   * Copy constructor.
   *
   * @param Cuckoo The source of information
   */
  cuckoo(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    copy(Cuckoo);
  }

  /**
   * Constructor from range [fisrt, last).
   * 
   * @param first The begin of range iterator
   * @param last The end of range iterator
   */
  template <class InputIterator>
  cuckoo(InputIterator first, InputIterator last, 
         const Hash& hasher = Hash(), 
         const Equal& equal = Equal()) {
    d_ = D;
    init_length_ = INIT_LENGTH;
    max_loop_ = MAX_LOOP;
    step_ = STEP;
    hasher_ = hasher;
    key_equal_ = equal;
    init();
    size_t last_pos = last.pos;
    for (iterator it = first; it.pos != last_pos; ++it) {
      add_new(*it);
    }
  }

  /**
   * Update parameter of cuckoo, deleting all the data from it.
   *
   * @warning For test only!!!
   */
  void set_up(size_t d = D, size_t init_length = INIT_LENGTH, 
              size_t max_loop = MAX_LOOP, double step = STEP) {
    clear_all();
    d_ = d;
    init_length_ = init_length;
    max_loop_ = max_loop;
    step_ = step;
    init();
  }

  /**
   * Operator==
   *
   * @param Cuckoo Object to be compared with
   */
  bool operator==(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    return (*this).data_ == Cuckoo.data_;
  }

  /**
   * @see operator==
   */
  bool operator!=(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    return !(*this == Cuckoo);
  }

  /**
   * Swap data with Cuckoo.
   *
   * @param Cuckoo Cuckoo of the same type
   */
  void swap(cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    std::swap(d_, Cuckoo.d_);
    std::swap(init_length_, Cuckoo.init_length_);
    std::swap(max_loop_, Cuckoo.max_loop_);
    std::swap(step_, Cuckoo.step_);
    std::swap(hasher_, Cuckoo.hasher_);
    std::swap(key_equal_, Cuckoo.key_equal_);
    std::swap(data_, Cuckoo.data_);
    std::swap(exists_, Cuckoo.exists_);
    std::swap(len_, Cuckoo.len_);
    std::swap(len_part_, Cuckoo.len_part_);
    std::swap(size_, Cuckoo.size_);
  }

  /**
   * Get iterator to begin of cuckoo.
   */
  inline iterator begin() {
    iterator it = iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }
  
  /**
   * Get const_iterator to begin of cuckoo (for const objects).
   */
  inline const_iterator begin() const {
    const_iterator it = const_iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }
  
  /**
   * Get iterator to end of cuckoo.
   */
  inline iterator end() {
    return iterator(len_, this);
  }

  /**
   * Get const_iterator to begin of cuckoo (for const objects).
   */
  inline const_iterator end() const {
    return const_iterator(len_, this);
  }

  /**
   * Get value by key or created pair key-value.
   *
   * @param k Key value
   */
  Value& operator[](const Key& k) {
    return (*((this->insert(make_pair(k, Value()))).first)).second;
  }

  /**
   * Erase data at iterator.
   *
   * @param it Iterator to Data element to be removed
   */
  void erase(iterator it) {
    remove(it);
  }

  /**
   * Erase range of Data elements.
   *
   * @param first The begin of range iterator
   * @param last The end of range iterator
   */
  void erase(iterator first, iterator last) {
    while (first.pos != last.pos) {
      first = remove(first);
    }
  }

  /**
   * Erase element by key.
   *
   * @param k Key value
   * @return Returns 1 if element was erased and 0 if it didn't exist
   */
  size_t erase(const Key& k) {
    iterator it = find(k);
    if (it.pos != len_) {
      remove(it);
      return 1;
    }
    return 0;
  }

  /**
   * Find element by key.
   *
   * @param k Key value
   * @return Returns iterator to element or to end of cuckoo, if element 
   * doesn't exist
   */
  iterator find(const Key& k) {
    size_t dist = 0;
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      size_t position = pos + dist;
      if (key_equal_(data_[i][pos].first, k) && (get_exists(position))) {    
        return iterator(position, this);
      }
      dist += len_part_;
    }
    return end();
  }
  
  /**
   * Find element by key (for const objects).
   *
   * @param k Key value
   * @return Returns const_iterator to element or to end of cuckoo, if element 
   * doesn't exist
   */
  const_iterator find(const Key& k) const {
    size_t dist = 0;
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      size_t position = pos + dist;
      if (key_equal_(data_[i][pos].first, k) && (get_exists(position))) {    
        return const_iterator(position, this);
      }
      dist += len_part_;
    }
    return end();
  }

  /**
   * Count number of elements with this key.
   *
   * @param k Key value
   * @return Returns 1 if element exists and 0 otherwise
   */
  size_t count(const Key& k) const {
    return (find(k)).pos != len_;
  }
  
  /**
   * Finds range of elements with key.
   *
   * @param k Key value
   * @return Returns pair that determines the range [fisrt, last) or
   * pair with both iterator pointing to the end of cuckoo.
   */ 
  pair<iterator, iterator> equal_range(const Key& k) {
    iterator l = find(k);
    iterator r = l;
    return std::make_pair<iterator, iterator>(l, ++r);
  }

  /**
   * Finds range of elements with key (for const objects).
   *
   * @param k Key value
   * @return Returns pair that determines the range [fisrt, last) or
   * pair with both iterator pointing to the end of cuckoo.
   */ 
  pair<const_iterator, const_iterator> equal_range(const Key& k) const {
    const_iterator l = find(k);
    const_iterator r = l;
    return std::make_pair<const_iterator, const_iterator>(l, ++r);
  }

  /** 
   * Inserts Data element to cuckoo.
   *
   * @param k The new Data element
   * @return Return pair with iterator to existing element and bool value,
   * which is true if element was inserted or false if it existed before.
   */
  pair<iterator, bool> insert(const Data& k) {
    iterator res = find(k.first);
    if (res.pos != len_) {
      return make_pair(res, false);
    } 
    add_new(k);
    return make_pair(iterator(hash(k.first, 0), this), true);
  }

  /**
   * Inserts range of Data elements.
   *
   * @param first The begin of range iterator
   * @param last The end of range iterator
   */
  template <class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    size_t last_pos = last.pos;
    for (iterator it = first; it.pos != last_pos; ++it) {
      insert(*it);
    }
  } 

  /** 
   * Clear all data from cuckoo.
   */
  void clear() {
    clear_all();
    init();
  }

  /**
   * Check whether cuckoo is empty.
   * 
   * @return Returns true of cuckoo is empty and false otherwise
   */
  inline bool empty() const {
    return (size_ == 0);
  }
 
  /**
   * Shows size of cuckoo.
   *
   * @return Returns size of cuckoo
   */ 
  inline size_t size() const {
    return size_;
  }

  /**
   * Shows length of cuckoo (actual number of elements).
   *
   * @return Returns length of cuckoo
   */ 
  inline size_t length() const {
    return len_;
  }
};

#endif /* _CUCKOO_HPP_ */

/*
 * trie.hpp
 *
 * Created on: 3.04.2011
 *     Author: Mariya Fomkina
 */

#ifndef _TRIE_HPP_
#define _TRIE_HPP_

using namespace std;

template <class Key, class Value>
class trie {
private:
  struct Node {
    // can be one of two values - 0 or 1 (binary alphabet)
    char letter;
    // these two are the numbers of cells with 0 and 1 next letters
    unsigned long zero; 
    unsigned long one;
    // the pointer to data stored
    Value* data; 
    Node() : letter(0), zero(0), one(0), data(NULL) {}
  };
  vector<Node> tree_;
  vector<Value> data_;
  size_t size_;
  // data storage should better be in arrays,
  // but it's possible to test on vectors
  //Node* tree_;
  //Value* data_;
public:
  class iterator {
  private: 
    size_t pos;
    trie* tree;
    iterator(size_t p, trie* t) : pos(p), tree(t) {}
    friend class trie;
  public:
    iterator() : pos(0), tree(NULL) {}

    void operator=(const iterator &it) {
      pos = it.pos;
      tree = it.tree;
    }

    iterator& operator++() {
      //some code
      return *this;
    }

    iterator operator++(int) {
      iterator res = *this;
      this->operator++();
      return res;
    }

    pair<Key, Value>& operator*() {
      
    }

    bool operator==(const iterator &it) {
      return pos == it.pos && tree == it.tree;
    }

    bool operator!=(const iterator &it) {
      return !(*this == it);
    }
  };

private:
  void init() {
    size_ = 0;
  }
 
public:

  trie() { 
    init();
  }

  ~trie() {
  
  }

  trie(trie<Key, Value>& Trie) {
    init();
    *this = Trie;
  }

  trie<Key, Value>& operator=(const trie<Key, Value>& Trie) {
    //some operations of copiing
    return *this;
  }

  Value& operator[](const Key& k) {
    iterator it = find(k);
    if (it == end()) {
      it = insert(make_pair(k, Value())).first;
    }
    return (*it).second;
  }

  iterator begin() {

  }

  iterator end() {

  }

  iterator find(const Key& k) {

  }

  pair<iterator, bool> insert(const pair<const Key, Value> &k) {

  }

  size_t erase(const Key& k) {

  }

  void clear() {

  }

  bool empty() {
    return (size_ == 0);
  }

  size_t size() const {
    return size_;
  }
};
#endif

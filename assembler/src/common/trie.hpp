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
    // these are the numbers of cells with 0 and 1 next letters
    size_t zero; 
    size_t one;
    // parent node
    size_t parent;
    // the pointer to data stored
    Value* data; 
    Node(char l) : letter(l), zero(0), one(0), data(NULL) {}
  };
  vector<Node> tree_;
  vector<Value> data_;
  size_t size_;
  size_t len_;
  // data storage should better be in arrays,
  // but it's possible to test on vectors
  //Node* tree_;
  //Value* data_;
public:
  class iterator {
  private: 
    size_t pos;
    trie* tree;
    Key key;
    iterator(size_t p, trie* t) : pos(p), tree(t) {}
    friend class trie;
  public:
    iterator() : pos(0), tree(NULL) {}

    void operator=(const iterator &it) {
      pos = it.pos;
      tree = it.tree;
    }

    iterator& operator++() {
      while (((*tree).tree_[pos].data == NULL) && (pos < *tree.len_)) {
        ++pos;
      }
      //update key - the most boring process here
      return *this;
    }

    iterator operator++(int) {
      iterator res = *this;
      this->operator++();
      return res;
    }

    pair<Key, Value> operator*() {
      return make_pair(key, *((*tree).tree_[pos].data));
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
    len_ = 1;
    tree_.push_back(Node(255));
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
    iterator it(0, *this);
    return ++it;
  }

  iterator end() {
    iterator it(len_, *this);
    return it;
  }

  iterator find(const Key& k) {
    //TEMP
    return end();
  }

  pair<iterator, bool> insert(const pair<const Key, Value> &k) {
    //TEMP
    return make_pair(end(), true); 
  }

  size_t erase(const Key& k) {
    //TEMP
    return 0;
  }

  void clear() {
    //TEMP
  }

  bool empty() {
    return (size_ == 0);
  }

  size_t size() const {
    return size_;
  }

  size_t len() const {
    return len_;
  }
};
#endif

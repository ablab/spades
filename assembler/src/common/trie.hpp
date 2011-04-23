/*
 * trie.hpp
 *
 * Created on: 3.04.2011
 *     Author: Mariya Fomkina
 *   Modified: 23.04.2011 by author
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
    // these are the numbers of cells with 0 and 1 next letters and
    // parent node
    size_t nodes[3];
    // the pointer to data stored
    Value* data; 
    Node(char l) : letter(l), data(NULL) {
      nodes[0] = 0;
      nodes[1] = 0;
      nodes[2] = 0;
    }
  };
  vector<Node> tree_;
  vector<Value> data_;
  size_t size_;
  size_t len_;
public:
  class iterator {
  private: 
    size_t pos;
    trie* tree;
    Key key;
    pair<Key, Value> cur_pair_;
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
      char* temp = new char[sizeof(Key) * 8];
      size_t p = pos;
      for (size_t i = sizeof(Key) * 8 - 1; i >= 0; --i) {
        temp[i] = tree_[pos].letter;
        p = tree_[pos].nodes[2];
      }
      char* key_temp = new char[sizeof(Key)];
      for (int i = 0; i < sizeof(Key); ++i) {
        key_temp[i] = temp[i * 8];
        for (size_t j = 1; j < 8; ++j) {
          key_temp[i] = (key_temp[i] << 1) + temp[i * 8 + j];
        }
      }
      memcpy(&key, key_temp, sizeof(Key));
      delete [] temp;
      delete [] key_temp;
      // the end of terrible block of code that 
      // should be rewritten
      return *this;
    }

    iterator operator++(int) {
      iterator res = *this;
      this->operator++();
      return res;
    }

    pair<Key, Value>& operator*() {
      cur_pair_ = make_pair(key, *((*tree).tree_[pos].data));
      return cur_pair_;
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
    tree_.push_back(Node(15));
  }
 
  void add_new(const pair<const Key, Value>& value) {
    Key k = value.first;
    size_t sok = sizeof(Key) * 8;
    char* temp_key = new char[sok];
    for (size_t i = 0; i < sok; ++i) {
      temp_key[i] = (k & (1 << (sok - 1 - i))) >> (sok - 1  -i);
    }
    size_t pos = 0;
    size_t i = -1;
    while (tree_[pos].nodes[(size_t)temp_key[i + 1]] != 0) {
      pos = tree_[pos].nodes[(size_t)temp_key[i + 1]];
      ++i;
    }
    for (size_t j = i + 1; j < sok; ++j) {
      Node node(temp_key[j]);
      tree_[pos].nodes[(size_t)temp_key[j]] = tree_.size();
      node.nodes[2] = pos;
      tree_.push_back(node);
      pos = tree_.size() - 1;
    }
    data_.push_back(value.second);
    tree_[pos].data = &(data_[data_.size() - 1]);
    ++size_;
    len_ = tree_.size();
  }

  iterator search_for(const Key& k) {
    size_t sok = sizeof(Key) * 8;
    char* temp_key = new char[sok];
    for (size_t i = 0; i < sok; ++i) {
      temp_key[i] = (k & (1 << (sok - 1 - i))) >> (sok - 1  -i);
    }
    size_t pos = 0;
    size_t i = -1;
    while (tree_[pos].nodes[(size_t)temp_key[i + 1]] != 0) {
      pos = tree_[pos].nodes[(size_t)temp_key[i + 1]];
      ++i;
    }
    if (i == sok - 1) {
      return iterator(pos, this);
    } else {
      return iterator(len_, this);
    }
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
    size_ = Trie.size_;
    len_ = Trie.len_;
    //not safe, can cause memory leak
    tree_.clear();
    tree_ = Trie.tree_;
    data_.clear();
    data_ = Trie.data_;
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
    iterator it(0, this);
    return ++it;
  }

  iterator end() {
    iterator it(len_, this);
    return it;
  }

  iterator find(const Key& k) {
    return search_for(k);
  }

  pair<iterator, bool> insert(const pair<const Key, Value> &k) {
    iterator it = search_for(k.first);
    if (it == end()) {
      add_new(k);
      return make_pair(end(), true);
    } else {
      return make_pair(it, false);
    } 
  }

  size_t erase(const Key& k) {
    iterator it = search_for(k);
    if (it != end()) {
      //simple and not optimal
      tree_[it.pos].data = NULL;
    }
    return 0;
  }

  void clear() {
    tree_.clear();
    data_.clear();
    size_ = 0;
    len_ = 1;
    tree_.push_back(Node(15));
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

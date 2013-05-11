#ifndef AHO_CORASICK_H_
#define AHO_CORASICK_H_

#include <map>
#include <string>
#include <vector>

#include "comparator.hpp"

class Node {
 public:
  Node(Node *fail_node = NULL) : fail(fail_node), output(NULL), word_index(-1) { }

  Node* getLink(char c) const;
  bool isTerminal() const;

  std::map<char, Node *> links;
  Node *fail;		//pointer to the end of the longest suffix
  Node *output;	//pointer to the terminal node recognized at this state
  int word_index;
};

typedef std::map<std::string *, std::vector<int>, Compare> seq2index_t;

class AhoCorasick {
 public:
  AhoCorasick() : root(new Node()) {};
  void cleanup();
  void addString(const std::string * str);
  void init();
  seq2index_t search(const std::string& str);
 private:
  const Node * go(const Node * current_state, char c);
  void isFound(const Node * current_state, int pos, seq2index_t& seq2index_match);
  void insert_match(std::string * seq, int pos, seq2index_t& seq2index_match);

  Node * root;
  std::vector<std::string *> patterns;
};

#endif /* AHO_CORASICK_H_ */

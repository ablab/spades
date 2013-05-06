#include "aho_corasick.hpp"

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <queue>

using namespace std;

Node* Node::getLink(char c) const {
  auto iter = links.find(c);
  return (iter != links.end()) ? iter->second : NULL;
}

bool Node::isTerminal() const {
  return word_index >= 0;
}

void AhoCorasick::addString(const string * str) {
  Node *current_node = root;
  for (size_t i = 0; i < str->length(); ++i) {
    Node *child_node = current_node->getLink((*str)[i]);
    if (!child_node) {
      child_node = new Node(root);
      current_node->links[(*str)[i]] = child_node;
    }
    current_node = child_node;
  }
  current_node->word_index = patterns.size();
  patterns.push_back(const_cast<std::string *>(str));
}

void AhoCorasick::cleanup() {
  std::queue<Node *> toDel;
  toDel.push(root);
  while (!toDel.empty()) {
    Node * current = toDel.front();
    toDel.pop();
    std::map<char, Node *> links = current->links;
    for (auto it = links.begin(); it != links.end(); ++it) {
      toDel.push(it->second);
    }
    delete current;
  }
}

//create fail (longest suffix) links in BFS manner
//and init output sets/links for each node
void AhoCorasick::init() {
  queue<Node *> q;
  q.push(root);
  while (!q.empty()) {
    Node *current_node = q.front();
    q.pop();
    for (auto iter = current_node->links.begin();
         iter != current_node->links.end(); ++iter) {
      const char symbol = iter->first;
      Node * const child = iter->second;
      q.push(child);

      Node *fail_node = current_node->fail;
      //iterate by suffix links until we find the suffix that is followed by desired symbol (1)
      //or reach the root which has NULL fail link (2)
      while (fail_node) {
        Node *fail_candidate = fail_node->getLink(symbol);
        if (fail_candidate) {
          child->fail = fail_candidate;
          break; // (1)
        }
        fail_node = fail_node->fail; //quit if NULL (2)
      }

      //http://www.cs.uku.fi/~kilpelai/BSA05/lectures/slides04.pdf
      //out(u):= out(u) U out(f(u));
      //This is done because the patterns recognized at f(u)(if
      //any), and only those, are proper suffixes of L(u), and shall
      //thus be recognized at state u also.
      //In my case I do not store a set, but must traverse .output links until NULL
      child->output = (child->fail->isTerminal()) ? child->fail : child->fail->output;
    }
  }
}

//iterate through fail links unless we (1) met the Node followed by
//desired char or 2) reached root
const Node * AhoCorasick::go(const Node * current_state, char c) {
  while (current_state) {
    Node *candidate = current_state->getLink(c);
    if (candidate) {
      current_state = candidate;
      return current_state; //(1)
    }
    current_state = current_state->fail;
  }
  return root; //(2)
}

void AhoCorasick::insert_match(std::string * seq, int pos) {
  std::map<std::string*, std::vector<int>, Compare>::iterator it;

  if (seq2index_match.end() != (it = seq2index_match.find(seq))) {
    it->second.push_back(pos);
  } else {
    std::vector<int> ind;
    ind.push_back(pos);
    seq2index_match.insert(make_pair(seq, ind));
  }
}

//check if node is terminal and
//traverse output links to find nay patterns
//that are recognized at this state too
void AhoCorasick::isFound(const Node * current_state, int pos) {

  if (current_state->isTerminal()) {
    insert_match(patterns[current_state->word_index], pos - patterns[current_state->word_index]->length() + 1);
  }

  const Node *terminal_node = current_state->output;
  //out(u):= out(u) U out(f(u)); => iterate through f(u), f(f(u)), etc
  //as if some pattern is recognized at f(u), it should be recognized at u too.
  while (terminal_node) {
    insert_match(patterns[terminal_node->word_index], pos - patterns[terminal_node->word_index]->length() + 1);
    terminal_node = terminal_node->output;
  }
}

void AhoCorasick::search(const std::string& str) {
  seq2index_match.clear();
  const Node * current_state = root;
  std::vector<int> res;
  for (int i = 0; i < (int) str.length(); ++i) {
    current_state = go(current_state, str[i]);
    isFound(current_state, i);
  }
}

std::map<std::string*, std::vector<int>, Compare> AhoCorasick::getMatch() const {
  return seq2index_match;
}


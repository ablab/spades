#pragma once

#define CONST_ALL(x) std::cbegin(x), std::cend(x)
#define ALL(x) std::begin(x), std::end(x)

#include <vector>
#include <string>
#include <queue>
#include <algorithm>

#undef NDEBUG
#include <cassert>
#include <unordered_map>
#include <iostream>
#include <unordered_set>
#include <numeric>
#include "utils.hpp"
#include <type_traits>
#include "utils/logger/log_writers.hpp"

class EmptyClass {};

template <typename T, typename U>
constexpr bool is_same_v = std::is_same<T, U>::value;

template <class GraphPointer, class Base = EmptyClass>
class AbstractGraphPointer : public Base {
 public:
  AbstractGraphPointer() = default;
  AbstractGraphPointer(const Base &gp) : Base{gp} {}
  AbstractGraphPointer(Base &&gp) : Base{gp} {}
  bool is_divergent() const {
    return crtp_this()->next().size() > 1;
  }

  bool is_convergent() const {
    return crtp_this()->prev().size() > 1;
  }

  std::vector<std::pair<GraphPointer, char>> next_pairs() const {
    std::vector<std::pair<GraphPointer, char>> result;
    for (const auto cur : crtp_this()->next()) {
      result.emplace_back(cur, cur.letter());
    }
    return result;
  }

  std::vector<std::pair<GraphPointer, char>> prev_pairs() const {
    std::vector<std::pair<GraphPointer, char>> result;
    for (const auto cur : crtp_this()->prev()) {
      result.emplace_back(cur, cur.letter());
    }
    return result;
  }

 private:
  GraphPointer *crtp_this() { return static_cast<GraphPointer *>(this); }
  const GraphPointer *crtp_this() const { return static_cast<const GraphPointer *>(this); }
};

template <class GraphPointer>
class ReversalGraphPointer : public AbstractGraphPointer<ReversalGraphPointer<GraphPointer>, GraphPointer> {
 public:
  using AbstractGraphPointer<ReversalGraphPointer<GraphPointer>, GraphPointer>::AbstractGraphPointer;
  std::vector<ReversalGraphPointer> next() const {
    auto result = GraphPointer::prev();
    return std::vector<ReversalGraphPointer>(CONST_ALL(result));
  }

  std::vector<ReversalGraphPointer> prev() const {
    auto result = GraphPointer::next();
    return std::vector<ReversalGraphPointer>(CONST_ALL(result));
  }
};

class Graph {
 public:
  class GraphPointer : public AbstractGraphPointer<Graph::GraphPointer> {
   public:
    GraphPointer() : edge_id_(-1), position_(-1), pgraph_{nullptr} {
      // INFO("Empty GP constructed");
      // empty pointer
    }

    // GraphPointer() = default;
    bool is_empty() const { return edge_id_ == size_t(-1); }

    GraphPointer(const Graph *pgraph, size_t edge_id, size_t position)
        : pgraph_{pgraph}, edge_id_{edge_id}, position_{position} {
      assert(edge_id_ < pgraph_->edges_.size());
      assert(position_ < pgraph_->edges_[edge_id_].size());
    }
    GraphPointer(const GraphPointer &) = default;
    GraphPointer(GraphPointer &&) = default;
    GraphPointer &operator=(const GraphPointer &) = default;
    GraphPointer &operator=(GraphPointer &&) = default;
    ~GraphPointer() noexcept = default;

    std::vector<size_t> edges() const {
      if (is_empty()) {
        return {size_t(-1)};
      } else {
        return {edge_id_};
      }
    }

    bool operator==(const GraphPointer &other) const {
      return (edge_id_ == other.edge_id_) && (position_ == other.position_);
    }

    char letter() const {
      assert(edge_id_ < pgraph_->edges_.size());
      assert(position_ < pgraph_->edges_[edge_id_].size());
      return pgraph_->edges_[edge_id_][position_];
    }

    bool is_divergent() const {
      return (position_ + 1 == pgraph_->edges_[edge_id_].size()) && (pgraph_->outgoing_[edge_id_].size() > 1);
    }

    bool is_convergent() const { return (position_ == 0) && (pgraph_->ingoing_[edge_id_].size() > 1); }

    std::vector<GraphPointer> prev() const {
      if (position_ == 0) {
        std::vector<GraphPointer> result;
        for (size_t i : pgraph_->ingoing_[edge_id_]) {
          result.emplace_back(pgraph_, i, pgraph_->edges_[i].size() - 1);
        }

        return result;
      } else {
        return {GraphPointer(pgraph_, edge_id_, position_ - 1)};
      }
    }

    std::vector<GraphPointer> next() const {
      if (position_ + 1 == pgraph_->edges_[edge_id_].size()) {
        std::vector<GraphPointer> result;
        for (size_t i : pgraph_->outgoing_[edge_id_]) {
          assert(pgraph_->edges_[i].size() > 0);  // We do not allow ""-edges
          result.emplace_back(pgraph_, i, 0);
        }

        return result;
      } else {
        return {GraphPointer(pgraph_, edge_id_, position_ + 1)};
      }
    }

   private:
    friend struct std::hash<Graph::GraphPointer>;
    friend std::ostream &operator<<(std::ostream &os, const Graph::GraphPointer &p);
    const Graph *pgraph_;
    size_t edge_id_;
    size_t position_;
  };

  GraphPointer get_pointer(size_t edge_id, size_t position) const { return GraphPointer(this, edge_id, position); }

  explicit Graph(const std::vector<std::string> &edges) : edges_{edges}, ingoing_(edges.size()), outgoing_(edges.size()) {

  }

  std::vector<GraphPointer> begins() const {
    std::vector<GraphPointer> result;
    for (size_t i = 0; i < edges_.size(); ++i) {
      result.push_back(get_pointer(i, 0));
    }
    return result;
  }

  explicit Graph(size_t k, const std::vector<std::string> &edges) {
    std::unordered_map<std::string, std::vector<size_t>> prefix2ids;
    std::unordered_map<std::string, std::vector<size_t>> suffix2ids;
    for (size_t i = 0; i < edges.size(); ++i) {
      const auto &edge = edges[i];
      assert(edge.size() > k);
      prefix2ids[edge.substr(0, k)].push_back(i);
      suffix2ids[edge.substr(edge.size() - k)].push_back(i);
    }

    for (std::string edge : edges) {
      std::string b = edge.substr(0, k);
      std::string e = edge.substr(edge.size() - k);
      if (prefix2ids.count(e)) {
        edge.resize(edge.size() - k);
      }
      edges_.push_back(std::move(edge));
      outgoing_.push_back(prefix2ids[e]);
      ingoing_.push_back(suffix2ids[b]);
    }

    assert(check_symmetry());
  }

  void collapse_empty_edge_(size_t id) {
    // After calling this method even once
    // Graph is not an actual graph anymore (it does not have proper vertices)
    auto remove_element = [](std::vector<auto> &v, const auto &value) {
      return v.erase(std::remove(v.begin(), v.end(), value), v.end());
    };

    auto append = [](std::vector<auto> &v, const auto &w) {
      return v.insert(v.end(), ALL(w));
    };

    assert(edges_[id] == "");
    for (size_t i : ingoing_[id]) {
      assert(i != id);  // Edge is not zero-length loop
      remove_element(outgoing_[i], id);
      append(outgoing_[i], outgoing_[id]);
    }

    for (size_t i : outgoing_[id]) {
      assert(i != id);  // Edge is not zero-length loop
      remove_element(ingoing_[i], id);
      append(ingoing_[i], ingoing_[id]);
    }

    outgoing_[id] = ingoing_[id] = {};
    assert(check_symmetry()); // TODO remove it after debug
  }

  size_t collapse_empty_edges_() {
    size_t count = 0;
    for (size_t i = 0; i < edges_.size(); ++i) {
      if (edges_[i] == "") {
        collapse_empty_edge_(i);
        ++count;
      }
    }
    return count;
  }


  struct Vertex {
    std::vector<size_t> ingoing, outgoing;
  };

  std::vector<Vertex> vertices() const {
    std::vector<char> outgoing_processed(edges_.size(), false);
    std::vector<Vertex> vertices;
    for (size_t i = 0; i < edges_.size(); ++i) {
      if (!outgoing_processed[i] && outgoing_[i].size() > 0) {
        const auto &outgoing = outgoing_[i];
        const auto &ingoing = ingoing_[outgoing[0]];
        vertices.push_back({ingoing, outgoing});
        for (size_t id : ingoing) {
          outgoing_processed[id] = true;
        }
      }
    }

    return vertices;
  }

  bool check_symmetry() const {
    auto is_in = [](const auto &value, const auto &container) {
      return std::find(CONST_ALL(container), value) != container.cend();
    };

    for (size_t i = 0; i < edges_.size(); ++i) {
      for (size_t j : outgoing_[i]) {
        if (!is_in(i, ingoing_[j])) {
          ERROR("ingoing[" << j << "] does not contain " << i);
          return false;
        }
      }
      for (size_t j : ingoing_[i]) {
        if (!is_in(i, outgoing_[j])) {
          ERROR("outgoing[" << j << "] does not contain " << i);
          return false;
        }
      }
    }

    return true;
  }

  bool check_vertices() const {
    for (const auto &vertex : vertices()) {
      for (size_t i : vertex.ingoing) {
        assert(outgoing_[i] == vertex.outgoing);
      }
      for (size_t i : vertex.outgoing) {
        assert(ingoing_[i] == vertex.ingoing);
      }
    }

    return true;
  }

  void collapse_equal_prefix_(const Vertex &vertex) {
    if (vertex.outgoing.size() < 2) {
      return;
    }
    std::unordered_map<char, std::vector<size_t>> char2edges;
    for (size_t i : vertex.outgoing) {
      char2edges[edges_[i][0]].push_back(i);
    }

    std::vector<size_t> outgoing;

    for (const auto &kv : char2edges) {
      const auto &ids = kv.second;
      size_t len = shared_prefix_len_(ids);
      if (ids.size() == 1) {
        outgoing.insert(outgoing.end(), CONST_ALL(ids));
        continue;
      }
      std::string new_edge = edges_[ids[0]].substr(0, len);
      size_t new_id = edges_.size();
      for (size_t i : ids) {
        edges_[i] = edges_[i].substr(len);
        ingoing_[i] = { new_id };
      }
      edges_.push_back(new_edge);
      outgoing_.push_back(ids);
      ingoing_.push_back(vertex.ingoing);
      outgoing.push_back(new_id);
    }

    for (size_t i : vertex.ingoing) {
      outgoing_[i] = outgoing;
    }
  }

  void collapse() {
    // Add empty edge in order to collapse left (ingoing from nowhere) tips
    std::vector<size_t> left_tips;
    for (size_t i = 0; i < edges_.size(); ++i) {
      if (ingoing_[i].size() == 0) {
        left_tips.push_back(i);
      }
    }
    size_t new_edge = edges_.size();
    edges_.push_back("");
    ingoing_.push_back({});
    outgoing_.push_back(left_tips);
    for (size_t i : left_tips) {
      ingoing_[i] = {new_edge};
    }
    assert(check_vertices());
    assert(check_symmetry());

    INFO(n_bases() << " bases before collapsing");
    for (size_t iter = 0; iter < 10; ++iter) {
      for (const auto &vertex : vertices()) {
        collapse_equal_prefix_(vertex);
      }
      assert(check_vertices());
      INFO(n_bases() << " bases after " << iter << " collapsing iteration");
    }

    size_t empty_edges = collapse_empty_edges_();
    INFO(empty_edges << " empty edges removed");
  }

  // TODO implement right tips collapsing

  size_t shared_prefix_len_(const std::vector<size_t> &ids) const {
    assert(ids.size() > 0);
    size_t L = edges_[ids[0]].size();
    for (size_t i : ids) {
      L = std::min(L, edges_[i].size());
    }

    for (size_t len = 0; len < L; ++len) {
      for (size_t i : ids) {
        if (edges_[i][len] != edges_[ids[0]][len]) {
          return len;
        }
      }
    }

    return L;
  }

  std::vector<GraphPointer> all() const {
    std::vector<GraphPointer> result;
    for (size_t i = 0; i < edges_.size(); ++i) {
      for (size_t pos = 0; pos < edges_[i].size(); ++pos) {
        result.push_back(get_pointer(i, pos));
      }
    }
    return result;
  }

  size_t n_edges() const { return edges_.size(); }

  size_t n_vertices() const {
    return vertices().size();
  }

  size_t n_bases() const {
    size_t count = 0;
    for (const auto &edge : edges_) {
      count += edge.size();
    }
    return count;
  }

 private:
  std::vector<std::string> edges_;
  std::vector<std::vector<size_t>> ingoing_;
  std::vector<std::vector<size_t>> outgoing_;
  friend struct std::hash<Graph::GraphPointer>;
};

class DBGraph {
public:
    class GraphPointer : public AbstractGraphPointer<DBGraph::GraphPointer> {
    public:
      GraphPointer() : edge_id_(-1), position_(-1), pgraph_{nullptr} {
        // INFO("Empty GP constructed");
        // empty pointer
      }

      // GraphPointer() = default;
      bool is_empty() const {
        return edge_id_ == size_t(-1);
      }

      void normalize_prefix_to_suffix_() {
        while (position_ < pgraph_->k_ && pgraph_->ingoing_[edge_id_].size()) {
          size_t old_id = edge_id_;
          size_t old_pos = position_;

          edge_id_ = pgraph_->ingoing_[edge_id_][0];
          position_ = pgraph_->edges_[edge_id_].size() - pgraph_->k_ + position_;

          // assert(pgraph_->edges_[old_id][old_pos] == pgraph_->edges_[edge_id_][position_]);
        }
      }

      void normalize_suffix_() {
        size_t old_id = edge_id_;
        size_t old_pos = position_;

        if (position_ >= pgraph_->edges_[edge_id_].size() - pgraph_->k_) {
            size_t new_edge_id = pgraph_->suffix_brothers_[edge_id_][0];
            position_ += (pgraph_->edges_[new_edge_id].size() - pgraph_->edges_[edge_id_].size());
            edge_id_ = new_edge_id;
        }
        //
        // assert(pgraph_->outgoing_[old_id] == pgraph_->outgoing_[edge_id_]);
        // for (size_t i = old_pos, j = position_; i < pgraph_->edges_[old_id].size() && j < pgraph_->edges_[edge_id_].size(); ++i, ++j) {
        //   assert(pgraph_->edges_[old_id][i] == pgraph_->edges_[edge_id_][j]);
        // }
      }

      GraphPointer(const DBGraph *pgraph, size_t edge_id, size_t position) : pgraph_{pgraph}, edge_id_{edge_id}, position_{position} {
        normalize_prefix_to_suffix_();
        // normalize_suffix_();  // Wwe should use only ONE kind of normalization
      }

        GraphPointer(const GraphPointer&) = default;
        GraphPointer(GraphPointer&&) = default;
        GraphPointer& operator=(const GraphPointer&) = default;
        GraphPointer& operator=(GraphPointer&&) = default;
        ~GraphPointer() noexcept = default;

        std::vector<size_t> edges() const {
          if (is_empty()) {
            return {size_t(-1)};
          }
          std::vector<size_t> result;
          result.push_back(edge_id_);

          const auto &edge = pgraph_->edges_[edge_id_];
          const size_t k = pgraph_->k_;
          if (position_ < k) {
            // FIXME use brothers arrays
            result.insert(result.end(), CONST_ALL(pgraph_->ingoing_[edge_id_]));
            // add brothers
            if (pgraph_->ingoing_[edge_id_].size()) {
              size_t in = pgraph_->ingoing_[edge_id_][0];
              result.insert(result.end(), CONST_ALL(pgraph_->outgoing_[in]));
            }
          }
          if (position_ >= edge.size() - k) {
            result.insert(result.end(), CONST_ALL(pgraph_->outgoing_[edge_id_]));
            // add brothers
            if (pgraph_->outgoing_[edge_id_].size()) {
              size_t out = pgraph_->outgoing_[edge_id_][0];
              result.insert(result.end(), CONST_ALL(pgraph_->ingoing_[out]));
            }
          }

          remove_duplicates(result);
          return result;
        }

        bool operator==(const GraphPointer &other) const {
            return (edge_id_ == other.edge_id_) && (position_ == other.position_);
        }

        char letter() const {
            return pgraph_->edges_[edge_id_][position_];
        }

        bool is_convergent() const {
          return prev().size() > 1;  // TODO implement it more efficiently
        }

        bool is_divergent() const {
            return position_ + 1 == pgraph_->edges_[edge_id_].size() && pgraph_->outgoing_[edge_id_].size() > 1;
            // // return next().size() > 1;
            // return true;
        }

        std::vector<GraphPointer> prev() const {
          if (position_ == 0) {
            assert(pgraph_->ingoing_[edge_id_].size() == 0);
            return {};
          } else if (position_ != pgraph_->k_ || pgraph_->ingoing_[edge_id_].size() == 0) {
            return {GraphPointer(pgraph_, edge_id_, position_ - 1)};
          } else {
            assert(position_ == pgraph_->k_);
            assert(pgraph_->ingoing_[edge_id_].size() > 0);
            std::vector<GraphPointer> result;
            for (size_t i : pgraph_->ingoing_[edge_id_]) {
              result.emplace_back(pgraph_, i,
                                  pgraph_->edges_[i].size() - pgraph_->k_ - 1);
            }

            std::unordered_set<char> letters;
            for (const auto &cur : result) {
              letters.insert(cur.letter());
            }
            if (letters.size() != result.size()) {
              ERROR(letters);
              ERROR(result);
              for (auto res : result) {
                ERROR(pgraph_->edges_[res.edge_id_]);
              }
            }
            assert(letters.size() == result.size());

            return result;
          }
        }

        std::vector<GraphPointer> next() const {
            if (position_ + 1 < pgraph_->edges_[edge_id_].size()) {
                // assert(position_ != pgraph_->k_ - 1 || pgraph_->prefix_brothers_[edge_id_].size() == 1);
                return { GraphPointer(pgraph_, edge_id_, position_ + 1)  };
            } else {
                // assert(position_ + 1 == pgraph_->edges_[edge_id_].size());
                std::vector<GraphPointer> result;
                for (size_t i : pgraph_->outgoing_[edge_id_]) {
                    result.emplace_back(pgraph_, i, pgraph_->k_);  // Vertices are k-mers
                }
                //
                // std::unordered_set<char> letters;
                // for (const auto &cur : result) {
                //     letters.insert(cur.letter());
                // }
                // if (letters.size() != result.size()) {
                //     ERROR(letters);
                //     ERROR(result);
                //     for (auto res : result) {
                //         ERROR(pgraph_->edges_[res.edge_id_]);
                //     }
                // }
                // assert(letters.size() == result.size());
                //
                return result;
            }
        }

    private:
        friend struct std::hash<DBGraph::GraphPointer>;
        friend std::ostream& operator<<(std::ostream &os, const DBGraph::GraphPointer &p);
        const DBGraph* pgraph_;
        size_t edge_id_;
        size_t position_;
    };

    Graph general_graph() const {
      return Graph(k_, edges_);
    }

    struct Path {
      GraphPointer begin;
      std::string turns;
      GraphPointer end;
    };

    Path trace_exact_sequence(const std::string &s) const {
      assert(s.length() > k_);

      // auto it = kp1mer_index_.find(s.substr(0, k_ + 1));
      // if (it == kp1mer_index_.cend()) {
      //   ERROR("First (k+1)-mer not found");
      //   return "";
      // }
      // kmer index is too large, use linear search instead;

      const auto kp1mer = s.substr(0, k_ + 1);
      auto begin = find_kp1mer(kp1mer);
      if (begin.is_empty()) {
        INFO("Initial k+1-mer not found!" << begin);
        return {GraphPointer(), "?", GraphPointer()};
      }

      auto cur = begin;
      assert(cur.letter() == s[0]);

      std::string turns;
      for (size_t i = 1; i < s.length(); ++i) {
        auto next_pairs = cur.next_pairs();
        bool flag = false;
        for (const auto &gp_l : next_pairs) {
          if (gp_l.second == s[i]) {
            if (next_pairs.size() > 1) {
              // divergence
              turns += gp_l.second;
            }
            cur = gp_l.first;
            flag = true;
            break;
          }
        }
        if (!flag) {
          ERROR(next_pairs);
          ERROR("Required letter: " << s[i]);
          ERROR(i << " " << s);
          return {begin, turns + "?", GraphPointer()};
        }
      }

      return {begin, turns, cur};
    }

    GraphPointer get_pointer(size_t edge_id, size_t position) const {
        return GraphPointer(this, edge_id, position);
    }

    GraphPointer find_kp1mer(const std::string &kp1mer) const {
      for (size_t i = 0; i < edges_.size(); ++i) {
        const auto &edge = edges_[i];
        size_t pos = edge.find(kp1mer);
        if (pos != edge.npos) {
          return get_pointer(i, pos);
        }
      }
      return GraphPointer();  // empty pointer
    }


    explicit DBGraph(size_t k, const std::vector<std::string> &edges) : k_{k}, edges_{edges} {
        std::unordered_map<std::string, std::vector<size_t>> prefix2ids;
        std::unordered_map<std::string, std::vector<size_t>> suffix2ids;
        for (size_t i = 0; i < edges_.size(); ++i) {
            const auto &edge = edges_[i];
            assert(edge.size() > k_);
            prefix2ids[edge.substr(0, k_)].push_back(i);
            suffix2ids[edge.substr(edge.size() - k_)].push_back(i);
        }

        outgoing_.resize(edges_.size());
        ingoing_.resize(edges_.size());
        for (size_t i = 0; i < edges_.size(); ++i) {
            const auto &edge = edges_[i];
            outgoing_[i] = prefix2ids[edge.substr(edge.size() - k_)];
            ingoing_[i] = suffix2ids[edge.substr(0, k_)];
            for (size_t idx : outgoing_[i]) {
                assert(edge.substr(edge.size() - k_) == edges_[idx].substr(0, k_));
            }
        }

        suffix_brothers_.resize(edges_.size());
        for (const auto &kv : suffix2ids) {
            assert(kv.first.size() == k_);
            const auto &indices = kv.second;  // Already sorted
            for (size_t i : indices) {
                suffix_brothers_[i] = indices;
            }
        }

        prefix_brothers_.resize(edges_.size());
        for (const auto &kv : prefix2ids) {
            assert(kv.first.size() == k_);
            const auto &indices = kv.second;  // Already sorted
            for (size_t i : indices) {
                prefix_brothers_[i] = indices;
            }
        }

        // // Construct (k+1)-mer index; (k+1)-mers must be unique!
        // for (size_t i = 0; i < edges_.size(); ++i) {
        //   const auto &edge = edges_[i];
        //   for (size_t j = 0; j < edge.size() - k_; ++j) {
        //     const auto kp1mer = edge.substr(j, j + k_ + 1);
        //     assert(kp1mer_index_.count(kp1mer) == 0);
        //     kp1mer_index_[kp1mer] = get_pointer(i, j);
        //   }
        // }
    }

    std::vector<GraphPointer> all() const;

    std::vector<size_t> walk(const std::vector<size_t> &edges, const size_t pathlen_, const bool forward) const {
        // for each edge we have pathlen to go (just a hastable map<edge, int>)
        // each edge has 3 states: unvisited, visited, and processed
        // we process all visited edges, for processing we mean visiting all outgoing (or ingoing) edges with pathlen = pathlen - current_edge_len + k - 1
        // all edges with updated pathlen become visited
        // we use heap to find visited (but not processed) edge with largest pathlen
        std::unordered_map<size_t, size_t> edge2pathlen;

        for (size_t edge_id : edges) {
            edge2pathlen[edge_id] = pathlen_;
        }

        std::priority_queue<std::pair<size_t, size_t>> q;
        // Fill the queue
        for (const auto &kv : edge2pathlen) {
            q.emplace(kv.second, kv.first); // TODO emplace pathlen - current_len instead
        }

        std::unordered_set<size_t> processed;

        const auto &further_edges = (forward) ? outgoing_ : ingoing_;

        while (!q.empty()) {
            auto kv = q.top(); q.pop();
            // Dirty heap idiom. We do not process one element twice untill it is not updated
            size_t edge_id = kv.second;
            if (!processed.count(edge_id)) {
                processed.insert(edge_id);  // Remove from the heap
            } else {
                continue;
            }

            size_t pathlen = edge2pathlen[edge_id];
            if (edges_[edge_id].size() < pathlen) {
                size_t extra_pathlen = pathlen - edges_[edge_id].size() + k_;
                for (size_t out_id : further_edges[edge_id]) {
                    if (edge2pathlen[out_id] < extra_pathlen) {
                        edge2pathlen[out_id] = extra_pathlen;
                        q.emplace(extra_pathlen, out_id);
                        processed.erase(out_id);  // Return to the heap
                    }
                }
            }
        }

        std::vector<size_t> result;
        for (const auto &kv : edge2pathlen) {
            if (kv.second > 0) {
                result.push_back(kv.first);
            }
        }

        remove_duplicates(result);

        return result;
    }


    DBGraph subgraph(const std::vector<size_t> &forward_edges,
                    const std::vector<size_t> &backward_edges,
                    size_t pathlen) {
        auto subgraph_forward = this->walk(forward_edges, pathlen, true);
        auto subgraph_backward = this->walk(backward_edges, pathlen, false);

        std::vector<size_t> subgraph_edges;
        subgraph_edges.insert(subgraph_edges.end(), CONST_ALL(subgraph_forward));
        subgraph_edges.insert(subgraph_edges.end(), CONST_ALL(subgraph_backward));
        remove_duplicates(subgraph_edges);

        std::vector<std::string> new_edges;
        new_edges.reserve(subgraph_edges.size());
        for (size_t i : subgraph_edges) {
            new_edges.push_back(edges_[i]);
        }

        return DBGraph(k_, new_edges);
    }


    size_t n_edges() const {
      return edges_.size();
    }

    size_t n_vertices() const {
      std::unordered_set<std::string> vertices;
      for (const auto &edge : edges_) {
        vertices.insert(edge.substr(0, k_));
        vertices.insert(edge.substr(edge.size() - k_));
      }
      return vertices.size();
    }

    // TODO implement ergo (condensation) loops detection
    bool loops() const {
      std::vector<char> visited(edges_.size(), 0); // 0 --- not visited, 1 --- current; 2 --- processed
      auto dfs = [this, &visited](size_t i, const auto& dfs) -> bool {
        // Pass lambda inside as a parameter. A trick from
        // https://stackoverflow.com/questions/2067988/recursive-lambda-functions-in-c11
        if (visited[i] == 2) {
          return false;
        }

        if (visited[i] == 1) {
          return true;
        }

        assert(visited[i] == 0);
        visited[i] = 1;
        for (size_t j : this->outgoing_[i]) {
          if (dfs(j, dfs)) {
            return true;
          }
        }

        visited[i] = 2;
        return false;
      };

      for (size_t i = 0; i < edges_.size(); ++i) {
        if (dfs(i, dfs)) {
          return true;
        }
      }
      return false;
    }

    DBGraph reversed() const {
      auto reversed_edges = edges_;
      for (auto &edge : reversed_edges) {
        std::reverse(ALL(edge));
      }
      return DBGraph(k_, reversed_edges);
    }


   private:
    size_t k_;
    std::vector<std::string> edges_;
    std::vector<std::vector<size_t>> ingoing_;
    std::vector<std::vector<size_t>> outgoing_;
    std::vector<std::vector<size_t>> suffix_brothers_;
    std::vector<std::vector<size_t>> prefix_brothers_;
    // std::unordered_map<std::string, GraphPointer> kp1mer_index_;
};


inline std::ostream& operator<<(std::ostream &os, const DBGraph::GraphPointer &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.edge_id_ << ", " << p.position_ << ")";  // << (p.empty() ? '@' : p.letter());
  }
}

inline std::ostream& operator<<(std::ostream &os, const Graph::GraphPointer &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.edge_id_ << ", " << p.position_ << ")";
  }
}

inline std::ostream& operator<<(std::ostream &os, const DBGraph::Path &path) {
    return os << "<" << path.begin << "-" << path.turns << "-" << path.end << ">";
}


inline size_t hash_size_t_pair(size_t s0, size_t s1) {
  s1 ^= s1 << 23; // a
  return ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) + s0;
}

namespace std {
template <>
struct hash<DBGraph::GraphPointer> {
  std::size_t operator()(const DBGraph::GraphPointer &p) const {
    return std::hash<size_t>()(hash_size_t_pair(p.edge_id_, p.position_));
  }
};

template <typename GraphPointer>
struct hash<ReversalGraphPointer<GraphPointer>> : public hash<GraphPointer> {};

template <>
struct hash<Graph::GraphPointer> {
  std::size_t operator()(const Graph::GraphPointer &p) const {
    return std::hash<size_t>()(hash_size_t_pair(p.edge_id_, p.position_));
  }
};
}  // namespace std

template <typename GraphPointer>
std::string restore_path(const GraphPointer &begin, const GraphPointer &end,
                         const std::vector<char> &turns) {
  if (begin.is_empty() || end.is_empty()) {
    return "";
  }

  auto cur = begin;

  size_t i = 0;
  std::string path;
  while (true) {
    path += cur.letter();

    if (cur == end && i == turns.size()) {
      break;
    }

    auto nexts = cur.next();
    assert(nexts.size());

    if (nexts.size() > 1) {
      cur = nexts[turns[i]];
      ++i;
    } else {
      cur = nexts[0];
    }
  }

  return path;
}

template <typename GraphPointer>
std::string restore_path(const GraphPointer &begin, const GraphPointer &end,
                         const std::string &turns) {
  if (begin.is_empty() || end.is_empty()) {
    return "";
  }

  auto cur = begin;

  size_t i = 0;
  std::string path;
  while (true) {
    path += cur.letter();

    if (cur == end && i == turns.size()) {
      break;
    }

    auto next_pairs = cur.next_pairs();
    assert(next_pairs.size());

    if (next_pairs.size() > 1) {
      bool flag = false;
      for (const auto &gp_l : next_pairs) {
        assert(i < turns.size());
        if (gp_l.second == turns[i]) {
          cur = gp_l.first;
          ++i;
          flag = true;
          break;
        }
      }
      if (!flag) {
        ERROR(path);
        ERROR("LETTER " << turns[i] << " not found");
        ERROR(next_pairs);
      }
      assert(flag);
    } else {
      cur = next_pairs[0].first;
    }
  }

  return path;
}

inline std::string restore_path(const DBGraph::Path &path) { return restore_path(path.begin, path.end, path.turns); }
// auto cur = graph.get_pointer(0, 0);

// std::cout << "Graph walk:" << std::endl;
// while (true) {
//     std::cout << cur.letter();
//     auto next = cur.next();
//     if (next.size()) {
//         cur = next[0];
//     } else {
//         break;
//     }
// }

#pragma once

#include "utils/logger/logger.hpp"

#include "llvm/ADT/IntrusiveRefCntPtr.h"

#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

enum class EventType { NONE, MATCH, INSERTION };
using score_t = double;


namespace pathtree {

struct Event {
  unsigned int m;  // We do not need size_t here
  EventType type;
};

template <class T>
class ObjectCounter {
 public:
  static size_t object_count_max() { return count_max_; }

  static size_t object_count_current() { return count_current_; }

  template <typename... Args>
  ObjectCounter(Args &&...) noexcept {
    ++count_current_;
    count_max_ = std::max(count_max_, count_current_);
  }

  ~ObjectCounter() { --count_current_; }

 private:
  static size_t count_max_;
  static size_t count_current_;
};

template <class T>
size_t ObjectCounter<T>::count_max_ = 0;

template <class T>
size_t ObjectCounter<T>::count_current_ = 0;

template <typename T>
class Node : public ObjectCounter<Node<T>>, public llvm::RefCountedBase<Node<T>> {
  using This = Node<T>;
  using ThisRef = llvm::IntrusiveRefCntPtr<This>;

 public:
  std::vector<T> collect() const {
    std::vector<T> result;
    const This *p = this;
    while (p) {
      result.push_back(p->payload_);
      p = p->parent_.get();
    }

    return result;
  }

  Node(const T &payload, const ThisRef &parent = nullptr) : payload_{payload}, parent_{parent} {}

  static ThisRef child(const T &payload, const ThisRef &parent = nullptr) { return new This(payload, parent); }

  const auto &payload() const { return payload_; }

 private:
  T payload_;
  ThisRef parent_;
};

template <class T>
using NodeRef = llvm::IntrusiveRefCntPtr<Node<T>>;

template <typename T>
NodeRef<T> make_child(const T &payload, const NodeRef<T> &parent = nullptr) {
  return Node<T>::child(payload, parent);
}

template <typename GraphCursor>
class PathLink : public llvm::RefCountedBase<PathLink<GraphCursor>> {
  using This = PathLink<GraphCursor>;
  using ThisRef = llvm::IntrusiveRefCntPtr<This>;

 public:
  double score() const {
    if (scores_.size() == 0) {
      return std::numeric_limits<double>::infinity();
    }

    return best_ancestor()->second.first;
  }

  auto best_ancestor() const {
    assert(scores_.size());
    return std::min_element(scores_.cbegin(), scores_.cend(),
                            [](const auto &e1, const auto &e2) { return e1.second.first < e2.second.first; });
  }

  bool update(GraphCursor gp, double score, const ThisRef &pl) {
    auto val = std::make_pair(score, pl);
    auto it_fl = scores_.insert({gp, val});
    bool inserted = it_fl.second;
    if (!inserted) {
      const auto &it = it_fl.first;
      if (it->second.first > score) {
        it->second = std::move(val);
        assert(scores_.size());
        return true;
      } else {
        assert(scores_.size());
        return false;
      }
    } else {
      assert(scores_.size());
      return true;
    }
  }

  static ThisRef master_source() {
    auto result(new This());
    result->scores_[GraphCursor()] = std::make_pair(0, nullptr);  // master_source score should be 0
    return result;
  }

  auto best_scores() const {
    struct Payload {
      const This *state;
      score_t score;
      const This *ancestor;
    };

    struct TopScore {
      score_t score;
      const This *ancestor;
    };

    std::unordered_map<const This *, TopScore> result;
    std::unordered_map<const This *, std::unordered_set<This *>> forward;

    struct Comp {
      bool operator()(const Payload &e1, const Payload &e2) const { return e1.score > e2.score; }
    };

    std::priority_queue<Payload, std::vector<Payload>, Comp> q;
    for (const auto &kv : scores_) {
      This *p = kv.second.second.get();
      q.push({p, kv.second.first, this});
      forward[p].insert(const_cast<This *>(this));
    }

    while (!q.empty()) {
      auto payload = q.top();
      q.pop();
      const This *state;
      score_t score;
      const This *ancestor;
      std::tie(state, score, ancestor) = {payload.state, payload.score, payload.ancestor};

      if (result.count(state)) {
        continue;
      }
      result[state] = {score, ancestor};

      auto best = state->best_ancestor();
      for (auto it = state->scores_.cbegin(); it != state->scores_.cend(); ++it) {
        score_t delta = it->second.first - best->second.first;
        This *p = it->second.second.get();
        if (p) {  // is not master source TODO make a method
          q.push({p, score + delta, state});
          forward[p].insert(const_cast<This *>(state));
        }
      }
    }

    return make_pair(result, forward);
  }

  auto exit_scores() const {
    std::unordered_map<const This *, score_t> result;
    for (const auto &kv : scores_) {
      result[kv.first] = kv.second.fisrt;
    }

    return result;
  }

  template <typename Function>
  auto apply(const Function &function) {
    std::unordered_set<const This *> checked;

    std::queue<const This *> q;
    q.push(this);
    while (!q.empty()) {
      const This *current = q.front();
      q.pop();
      if (!current || checked.count(current)) {
        continue;
      }

      function(*current);

      checked.insert(current);

      for (const auto &kv : current->scores_) {
        const This *p = kv.second.second.get();
        if (p) {
          q.push(p);
        }
      }
    }
  }

  auto states_without_children() const {
    std::unordered_set<const This *> checked;
    std::unordered_set<const This *> result;

    std::queue<const This *> q;
    q.push(this);
    while (!q.empty()) {
      const This *current = q.front();
      q.pop();
      if (!current || checked.count(current)) {
        continue;
      }

      if (current->scores_.size() == 0) {
        result.insert(current);
      }

      checked.insert(current);

      for (const auto &kv : current->scores_) {
        const This *p = kv.second.second.get();
        if (p) {
          q.push(p);
        }
      }
    }

    return result;
  }

  bool check_all_states_have_children() const { return states_without_children().size() == 0; }

  void clean_non_aggressive() {
    if (!check_all_states_have_children()) {
      auto bad_states = states_without_children();
      INFO(this);
      INFO(bad_states);
      auto fnc = [&](const auto &state) {
        for (const auto &kv : state.scores_) {
          const This *p = kv.second.second.get();
          if (bad_states.count(p)) {
            ERROR(kv.first.is_empty());
            ERROR(kv.second.first);
            ERROR(kv.first.letter());
          }
        }
      };

      apply(fnc);
    }
    assert(check_all_states_have_children());

    auto bs_fwd = best_scores();
    auto bs = bs_fwd.first;
    auto forward = bs_fwd.second;

    std::queue<This *> q;
    std::unordered_set<This *> processed;
    std::unordered_set<This *> collapsed;
    q.push(this);

    while (!q.empty()) {
      This *current = q.front();
      q.pop();

      if (processed.count(current)) {
        continue;
      }
      processed.insert(current);

      // The order of the elements that are not erased is preserved (this makes it possible to erase individual elements
      // while iterating through the container)
      for (auto it = current->scores_.begin(); it != current->scores_.end();) {
        const This *p = it->second.second.get();
        assert(p);
        assert(bs.count(p));

        if ((current == this && bs[p].ancestor != this) || (current != this && bs[p].ancestor == this)) {
          it = current->scores_.erase(it);
          TRACE("Link erased");
        } else {
          ++it;
        }
      }

      if (current->scores_.size() == 0) {
        TRACE("We collapse current state");
        collapsed.insert(current);
        continue;
      }

      current->clean_left_link_non_aggressive();

      for (auto it = current->scores_.begin(); it != current->scores_.end(); ++it) {
        This *p = it->second.second.get();
        // We cannot do this in the previous loop since some references could be invalidated
        if (!it->first.is_empty()) {
          q.push(p);
        }
      }
    }

    assert(q.empty());
    // Remove collapsed states recursively
    for (This *p : collapsed) {
      for (This *pp : forward[p]) {
        q.push(pp);
      }
    }

    while (!q.empty()) {
      This *current = q.front();
      q.pop();

      for (auto it = current->scores_.begin(); it != current->scores_.end();) {
        This *p = it->second.second.get();
        assert(p);

        if (collapsed.count(p)) {
          it = current->scores_.erase(it);
          TRACE("Link erased");
        } else {
          ++it;
        }
      }

      if (current->scores_.size() == 0) {
        TRACE("We collapse current state");
        collapsed.insert(current);
        for (This *pp : forward[current]) {
          q.push(pp);
        }
      }
    }
    assert(q.empty());
    assert(check_all_states_have_children());
    assert(!collapsed.count(this));

    INFO(collapsed.size() << " states collapsed");
  }

  struct AnnotatedPath {
    std::vector<GraphCursor> path;
    double score;
    std::vector<Event> events;
  };

  std::vector<AnnotatedPath> top_k(size_t k) const {
    // using Payload = std::tuple<GraphCursor, double, const This *>;
    struct Payload {
      GraphCursor gp;
      double cost;
      const This *p;
    };
    using Node = Node<Payload>;
    using SPT = NodeRef<Payload>;
    struct Comp {
      bool operator()(const SPT &e1, const SPT &e2) const {
        return e1->payload().cost> e2->payload().cost;
      }
    };
    std::priority_queue<SPT, std::vector<SPT>, Comp> q;
    auto extract_path = [&q]() -> AnnotatedPath {
      assert(!q.empty());
      SPT tail = q.top();
      q.pop();
      // GraphCursor gp;
      const This *p = tail->payload().p;
      double cost = tail->payload().cost;

      TRACE("Extracting path with cost " << cost);

      while (p) {
        auto new_tail = tail;
        auto best = p->best_ancestor();
        for (auto it = p->scores_.cbegin(); it != p->scores_.cend(); ++it) {
          double delta = it->second.first - best->second.first;
          auto spt = make_child<Payload>({it->first, cost + delta, const_cast<const This *>(it->second.second.get())}, tail);
          assert(spt->payload().gp.is_empty() == (spt->payload().p->event.type == EventType::NONE));
          if (it != best) {
            q.push(spt);
          } else {
            new_tail = spt;
          }
        }
        p = best->second.second.get();
        tail = new_tail;

        if (p->best_ancestor()->second.second == nullptr) {  // TODO Add method to check master_source
          break;
        }
      }

      std::vector<GraphCursor> path;
      std::vector<Event> events;

      for (const auto &tpl : tail->collect()) {
        path.push_back(tpl.gp);
        events.push_back(tpl.p->event);
      }

      return AnnotatedPath{path, cost, events};
    };

    std::vector<AnnotatedPath> result;
    auto best = best_ancestor();
    if (best == scores_.end()) {
      return result;
      // TODO Support empty Link as a common case and remove this workaround
    }

    double best_score = best->second.first;
    auto initial = make_child<Payload>({GraphCursor(), best_score, this});
    q.push(initial);

    for (size_t i = 0; i < k && !q.empty(); ++i) {
      result.push_back(extract_path());
      TRACE((i + 1) << " top paths extracted");
      TRACE(Node::object_count_current() << " current # of best path tree nodes");
      TRACE(Node::object_count_max() << " max # of best path tree nodes");
    }

    return result;
  }

  void clean_left_link_non_aggressive() {
    auto it_terminal = scores_.find(GraphCursor());
    if (it_terminal == scores_.end()) {
      return;
    }

    auto best = best_ancestor();
    if (best == scores_.end()) {
      return;
      // TODO Support empty Link as a common case and remove this workaround
    }

    if (best == it_terminal) {
      for (auto it = scores_.begin(); it != scores_.end();) {
        if (!it->first.is_empty()) {
          it = scores_.erase(it);
        } else {
          ++it;
        }
      }
    } else {
      scores_.erase(it_terminal);
    }
  }

  size_t clean_left_link_old_() {
    // the idea is in the following:
    // if a path is a prefix or suffix of another one,
    // we do not consider them as DIFFERENT paths. We leave only one (the best) of them
    // In order to implement this we do not allow to link to be "left-terminal" (contain backref to "empty cursor") and
    // "transit" (contain backrefs to other nontrivial links).
    // Here we just remove ref to empty cursor or refs to nontrivial links dependent on what is better

    // If terminal ref is present:
    // 1) Remove all other refs worse than it
    // 2) In case of some non-terminal refs still present, remove terminal ref as well
    auto it_terminal = scores_.find(GraphCursor());
    if (it_terminal == scores_.end()) {
      return 0;
    }

    double terminal_score = it_terminal->second.first;
    std::vector<GraphCursor> to_remove;
    for (const auto &kv : scores_) {
      if (kv.first.is_empty()) {
        continue;
      }
      double score = kv.second.first;
      if (score >= terminal_score) {  // We prefer short paths
        to_remove.push_back(kv.first);
      } else {
        to_remove.push_back(GraphCursor());
      }
    }

    size_t count = 0;
    for (const auto &key : to_remove) {
      count += scores_.erase(key);
    }

    return count;
  }

  size_t clean_left_link_() {
    // the idea is in the following:
    // if a path is a prefix or suffix of another one,
    // we do not consider them as DIFFERENT paths. We leave only one (the best) of them
    // In order to implement this we do not allow to link to be "left-terminal" (contain backref to "empty cursor") and
    // "transit" (contain backrefs to other nontrivial links).
    // Here we just remove ref to empty cursor or refs to nontrivial links dependent on what is better

    // If terminal ref is present:
    // 1) Remove all other refs worse than it
    // 2) In case of some non-terminal refs still present, remove terminal ref as well
    // We add some tolerance. If two paths have almost equal scores we keep them both
    const double eps = 1e-7;
    auto it_terminal = scores_.find(GraphCursor());
    if (it_terminal == scores_.end()) {
      return 0;
    }

    double terminal_score = it_terminal->second.first;
    std::vector<GraphCursor> to_remove;
    for (const auto &kv : scores_) {
      if (kv.first.is_empty()) {
        continue;
      }
      double score = kv.second.first;
      if (score > terminal_score + eps) {
        to_remove.push_back(kv.first);
      } else if (score < terminal_score - eps) {
        to_remove.push_back(GraphCursor());
      }
    }

    size_t count = 0;
    for (const auto &key : to_remove) {
      count += scores_.erase(key);
    }

    return count;
  }

  ThisRef clone() const { return new This(*this); }

  Event event;
 private:
  std::unordered_map<GraphCursor, std::pair<double, ThisRef>> scores_;
};

template <class GraphCursor>
using PathLinkRef = llvm::IntrusiveRefCntPtr<PathLink<GraphCursor>>;
}  // namespace pathtree

template <class GraphCursor>
class PathSet {
 public:
  using AnnotatedPath = typename pathtree::PathLink<GraphCursor>::AnnotatedPath;
  class path_container {
   public:
    path_container(const pathtree::PathLink<GraphCursor> &paths, size_t k) : paths_(paths.top_k(k)) {}

    auto begin() const { return paths_.begin(); }
    auto end() const { return paths_.end(); }
    size_t size() const { return paths_.size(); }
    bool empty() const { return paths_.empty(); }
    auto operator[](size_t n) const { return paths_[n]; }

    static std::string str(const std::vector<GraphCursor> &path) {
      std::string s;
      for (size_t i = 0; i < path.size(); ++i) {
        if (path[i].is_empty()) continue;

        s += path[i].letter();
      }
      return s;
    }

    static std::string compress_alignment(const std::string &s) {
      size_t count = 0;
      char prev_c = '\0';
      std::string result;
      for (size_t i = 0; i <= s.size(); ++i) {
        char c = s[i];
        if (c == prev_c) {
          ++count;
        } else {
          if (prev_c != '\0') {
            if (count > 1) {
              result += std::to_string(count);
            }
            if (prev_c == '-') {
              prev_c = 'D';
            }
            result += prev_c;
          }
          count = 1;
        }
        prev_c = c;
      }

      return result;
    }

    static std::string alignment(const AnnotatedPath &apath, size_t m) {
      std::string s;
      size_t prev_position = 0;
      assert(apath.path.size() == apath.events.size());
      for (size_t i = 0; i < apath.path.size(); ++i) {
        if (apath.path[i].is_empty()) {
          // assert(apath.events[i].type == EventType::NONE);
          continue;
        };

        if (apath.events[i].type == EventType::NONE) {
          ERROR("Position: " << i);
          ERROR("Path: " << str(apath.path));
          ERROR("Alignment:" << s);
        }
        assert(apath.events[i].type != EventType::NONE);
        for (size_t j = prev_position + 1; j < apath.events[i].m; ++j) {
          s += '-';
        }
        prev_position = apath.events[i].m;
        s += apath.events[i].type == EventType::MATCH ? 'M' : 'I';
      }

      // Add trailing gaps (-)
      for (size_t j = prev_position + 1; j <= m; ++j) {
          s += '-';
      }
      return s;
    }

    std::string str(size_t n) const { return str(paths_[n].path); }
    std::string alignment(size_t n, size_t m = 0) const { return alignment(paths_[n], m); }

   private:
    std::vector<AnnotatedPath> paths_;
  };

  pathtree::PathLink<GraphCursor> &pathlink() { return pathlink_; }
  const pathtree::PathLink<GraphCursor> &pathlink() const { return pathlink_; }

  double best_score() const { return pathlink_.score(); }
  AnnotatedPath best_path() const { return top_k(1)[0]; }
  std::string best_path_string() const { return path_container::str(best_path().path); }

  path_container top_k(size_t k) const { return path_container(pathlink_, k); }

  auto clip_tails_non_aggressive() { return pathlink_.clean_non_aggressive(); }

 private:
  pathtree::PathLink<GraphCursor> pathlink_;
};

// vim: set ts=2 sw=2 et :

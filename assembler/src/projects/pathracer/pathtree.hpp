#pragma once

#include "pathtrie.hpp"
#include "object_counter.hpp"

#include "utils/logger/logger.hpp"

#include <llvm/ADT/IntrusiveRefCntPtr.h>
#include <debug_assert/debug_assert.hpp>

#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fees.hpp"

enum class EventType { NONE, MATCH, INSERTION };
using score_t = double;

struct pathtree_assert : debug_assert::default_handler,
                         debug_assert::set_level<1> {};

namespace pathtree {

// TODO Pack Event (shoud weigth 4 bytes only!)
struct Event {
  unsigned int m;  // We do not need size_t here
  EventType type;
};

template <typename GraphCursor>
class PathLink : public llvm::RefCountedBase<PathLink<GraphCursor>>,
                 public AtomicObjectCounter<PathLink<GraphCursor>> {
  using This = PathLink<GraphCursor>;
  using ThisRef = llvm::IntrusiveRefCntPtr<This>;

  void* operator new (size_t sz) {
      return ::operator new(sz);
  }

 public:
  auto best_ancestor() const {
    DEBUG_ASSERT(scores_.size(), pathtree_assert{});
    return std::min_element(scores_.cbegin(), scores_.cend(),
                            [](const auto &e1, const auto &e2) { return e1.second.first < e2.second.first; });
  }

  double score() const {
    if (scores_.size() == 0) {
      return std::numeric_limits<double>::infinity();
    }

    return best_ancestor()->second.first;
  }

  bool update(GraphCursor gp, double score, const ThisRef &pl) {
    auto val = std::make_pair(score, pl);
    auto it_fl = scores_.insert({gp, val});
    bool inserted = it_fl.second;
    if (!inserted) {
      const auto &it = it_fl.first;
      if (it->second.first > score) {
        it->second = std::move(val);
        DEBUG_ASSERT(scores_.size(), pathtree_assert{});
        return true;
      } else {
        DEBUG_ASSERT(scores_.size(), pathtree_assert{});
        return false;
      }
    } else {
      DEBUG_ASSERT(scores_.size(), pathtree_assert{});
      return true;
    }
  }

  static ThisRef create() { return new This(); }
  ThisRef clone() const { return new This(*this); }

  static ThisRef master_source() {
    ThisRef result = create();
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
      std::tie(state, score, ancestor) = std::make_tuple(payload.state, payload.score, payload.ancestor);

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

  void clean_non_aggressive(const void *context) {
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
            ERROR(kv.first.letter(context));
          }
        }
      };

      apply(fnc);
    }
    DEBUG_ASSERT(check_all_states_have_children(), pathtree_assert{}, debug_assert::level<2>{});

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
        DEBUG_ASSERT(p, pathtree_assert{});
        DEBUG_ASSERT(bs.count(p), pathtree_assert{}, debug_assert::level<2>{});

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

    DEBUG_ASSERT(q.empty(), pathtree_assert{});
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
        DEBUG_ASSERT(p, pathtree_assert{});

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
    DEBUG_ASSERT(q.empty(), pathtree_assert{});
    DEBUG_ASSERT(check_all_states_have_children(), pathtree_assert{}, debug_assert::level<2>{});
    DEBUG_ASSERT(!collapsed.count(this), pathtree_assert{});

    INFO(collapsed.size() << " states collapsed");
  }

  struct AnnotatedPath {
    std::vector<GraphCursor> path;
    double score;
    std::vector<Event> events;

    bool empty() const { return path.empty(); }

    size_t size() const {
      DEBUG_ASSERT(path.size() == events.size(), pathtree_assert{});
      return path.size();
    }
  };

  std::vector<AnnotatedPath> top_k(size_t k) const {
    struct Event {
      GraphCursor gp; // PathLink does not know its own position!
      const This *path_link;
    };

    using EventPath = NodeRef<Event>;
    struct QueueElement {
      EventPath path;
      double cost;  // TODO Use scores as in the paper
    };

    struct Comp {
      bool operator()(const QueueElement &e1, const QueueElement &e2) const {
        return e1.cost > e2.cost;
      }
    };

    std::priority_queue<QueueElement, std::vector<QueueElement>, Comp> q;

    std::vector<AnnotatedPath> result;
    auto best = best_ancestor();
    if (best == scores_.end()) {
      return result;
      // TODO Support empty Link as a common case and remove this workaround
    }

    double best_score = best->second.first;
    auto SinkPath = make_root<Event>({GraphCursor(), this});
    q.push({SinkPath, best_score});

    auto get_annotated_path = [&](const EventPath &epath, double cost) -> AnnotatedPath {
      std::vector<GraphCursor> path;
      std::vector<pathtree::Event> events;

      for (const auto &tpl : epath->collect()) {
        if (tpl.gp.is_empty()) {
          continue;  // TODO Think about a better way to exclude empty cursors
        }
        path.push_back(tpl.gp);
        events.push_back(tpl.path_link->event);
      }

      return AnnotatedPath{path, cost, events};
    };

    while (!q.empty() && result.size() < k) {
      auto qe = q.top();
      q.pop();
      const GraphCursor &gp = qe.path->payload().gp;
      const This *path_link = qe.path->payload().path_link;
      const double &cost = qe.cost;

      TRACE("Extracting path with cost " << cost);
      // Check if path started with SOURCE: TODO implement it properly
      if (gp.is_empty() && !qe.path->is_root()) {
        // Report path
        auto annotated_path = get_annotated_path(qe.path, qe.cost);
        if (annotated_path.empty()) {
          // We should stop after the first empty path found
          // FIXME is it possible to obtain an empty path at all?????
          break;
        }
        result.push_back(get_annotated_path(qe.path, qe.cost));
        continue;
      }

      auto best = path_link->best_ancestor();
      for (auto it = path_link->scores_.cbegin(); it != path_link->scores_.cend(); ++it) {
        double delta = it->second.first - best->second.first;
        Event new_event{it->first, const_cast<const This *>(it->second.second.get())};
        auto new_path = qe.path->add(new_event);
        q.push({new_path, cost + delta});
      }
    }

    return result;
  }

  std::vector<AnnotatedPath> top_k_old(size_t k) const {
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
      DEBUG_ASSERT(!q.empty(), pathtree_assert{});
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
          DEBUG_ASSERT(spt->payload().gp.is_empty() == (spt->payload().p->event.type == EventType::NONE), pathtree_assert{});
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
        if (tpl.gp.is_empty()) {
          continue;  // TODO Think about a better way to exclude empty cursors
        }
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

    for (size_t i = 0; i < k && !q.empty();) {
      auto annotated_path = extract_path();
      if (annotated_path.empty()) {
        continue;
      }
      result.push_back(annotated_path);
      ++i;
      TRACE(i << " top paths extracted");
      TRACE(Node::object_count_current() << " current # of best path tree nodes");
      TRACE(Node::object_count_max() << " max # of best path tree nodes");
      TRACE(Node::object_count_constructed() << " overall best path tree nodes constructed");
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

    template <class Cursor>
    static std::string str(const std::vector<Cursor> &path,
                           const void *context) {
      std::string s;
      for (size_t i = 0; i < path.size(); ++i) {
        DEBUG_ASSERT(!path[i].is_empty(), pathtree_assert{});

        s += path[i].letter(context);
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
            result += std::to_string(count);
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

    static std::string alignment(const AnnotatedPath &apath, const hmm::Fees &fees,
                                 const void *context) {
      size_t m = fees.M;
      std::string s;
      size_t prev_position = 0;
      DEBUG_ASSERT(apath.path.size() == apath.events.size(), pathtree_assert{});
      for (size_t i = 0; i < apath.path.size(); ++i) {
        DEBUG_ASSERT(!apath.path[i].is_empty(), pathtree_assert{});

        if (apath.events[i].type == EventType::NONE) {
          ERROR("Position: " << i);
          ERROR("Path: " << str(apath.path, context));
          ERROR("Alignment:" << s);
        }
        DEBUG_ASSERT(apath.events[i].type != EventType::NONE, pathtree_assert{});
        for (size_t j = prev_position + 1; j < apath.events[i].m; ++j) {
          s += '-';
        }
        prev_position = apath.events[i].m;
        s += apath.events[i].type == EventType::MATCH ? (fees.consensus[apath.events[i].m - 1] == apath.path[i].letter(context) ? 'M' : 'X') : 'I';
      }

      // Add trailing gaps (-)
      for (size_t j = prev_position + 1; j <= m; ++j) {
          s += '-';
      }
      return s;
    }

    std::string str(size_t n, const void *context) const { return str(paths_[n].path, context); }
    std::string alignment(size_t n, const hmm::Fees &fees, const void *context) const { return alignment(paths_[n], fees, context); }

   private:
    std::vector<AnnotatedPath> paths_;
  };

  pathtree::PathLink<GraphCursor> &pathlink() { return pathlink_; }
  const pathtree::PathLink<GraphCursor> &pathlink() const { return pathlink_; }

  double best_score() const { return pathlink_.score(); }
  AnnotatedPath best_path() const { return top_k(1)[0]; }
  std::string best_path_string() const { return path_container::str(best_path().path); }

  path_container top_k(size_t k) const { return path_container(pathlink_, k); }

  auto clip_tails_non_aggressive(const void *context) { return pathlink_.clean_non_aggressive(context); }

 private:
  pathtree::PathLink<GraphCursor> pathlink_;
};

// vim: set ts=2 sw=2 et :

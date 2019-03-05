#pragma once

#include "cursor.hpp"
#include "fees.hpp"
#include "pathtree.hpp"
#include "depth_filter.hpp"

#include "utils/logger/logger.hpp"

#include <llvm/ADT/iterator.h>
#include <llvm/ADT/iterator_range.h>
#include <debug_assert/debug_assert.hpp>

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <memory>
#include <vector>
#include <limits>


extern "C" {
#include "hmmer.h"
}

namespace impl {

struct hmmpath_assert : debug_assert::default_handler,
                        debug_assert::set_level<1> {};

using pathtree::PathLink;
using pathtree::PathLinkRef;

template <typename Map>
class FilterMapMixin {
 public:
  template <typename Predicate>
  size_t filter_key(const Predicate &predicate) {
    size_t count = 0;
    for (auto it = crtp_this()->begin(); it != crtp_this()->end();) {
      if (predicate(it->first)) {
        it = crtp_this()->erase(it);
        ++count;
      } else {
        ++it;
      }
    }

    return count;
  }

  template <typename Predicate>
  size_t filter_key_value(const Predicate &predicate) {
    size_t count = 0;
    for (auto it = crtp_this()->begin(); it != crtp_this()->end();) {
      if (predicate(*it)) {
        it = crtp_this()->erase(it);
        ++count;
      } else {
        ++it;
      }
    }

    return count;
  }

 private:
  Map *crtp_this() { return static_cast<Map *>(this); }
  const Map *crtp_this() const { return static_cast<const Map *>(this); }
};

template <typename Map>
class ScoresFilterMapMixin : public FilterMapMixin<Map> {
 public:
  std::vector<score_t> scores() const {
    std::vector<score_t> scores;
    scores.reserve(crtp_this()->size());
    for (auto &kv : *crtp_this()) {
      scores.push_back(Map::score_fnc(kv));
    }

    return scores;
  }

  size_t score_filter(size_t n, score_t score) {
    n = std::min(n, crtp_this()->size());
    if (n == 0) {
      crtp_this()->clear();
      return 0;
    }

    {
      auto scores = crtp_this()->scores();
      std::nth_element(scores.begin(), scores.begin() + n - 1, scores.end());
      score = std::min(score, scores[n - 1]);
    }

    auto pred = [score](const auto &kv) -> bool {
      return Map::score_fnc(kv) > score;
    };

    return crtp_this()->filter_key_value(pred);
  }

 private:
  Map *crtp_this() { return static_cast<Map *>(this); }
  const Map *crtp_this() const { return static_cast<const Map *>(this); }
};


template <typename GraphCursor>
struct ScoredPLink {
  PathLinkRef<GraphCursor> plink;
  score_t score;
};

template <typename GraphCursor>
struct State {
  const GraphCursor &cursor;
  const PathLinkRef<GraphCursor> &plink;
  score_t score;
};

template <typename Map, typename Iterator>
class StateIterator : public llvm::iterator_facade_base<StateIterator<Map, Iterator>,
                                                        std::forward_iterator_tag,
                                                        State<typename Map::key_type>> {
  using This = StateIterator<Map, Iterator>;
  using GraphCursor = typename Map::key_type;

 public:
  StateIterator(const Iterator &it) : it_{it} {}

  State<GraphCursor> operator*() const { return Map::kv2state(*it_); }

  This &operator++() {
    ++it_;
    return *this;
  }

  bool operator==(const This &that) const { return it_ == that.it_; }

 private:
  Iterator it_;
};

template <typename Map, typename Iterator>
auto make_state_iterator(const Map &, const Iterator &it) {
  return StateIterator<Map, Iterator>(it);
}

template <typename Map, typename Iterator>
class Key2KeyValueIterator : public llvm::iterator_facade_base<Key2KeyValueIterator<Map, Iterator>,
                                                               std::forward_iterator_tag,
                                                               typename Map::value_type> {
  using This = Key2KeyValueIterator<Map, Iterator>;
 public:
  Key2KeyValueIterator(const Map &map, const Iterator &it) : map_{map}, it_{it} {}

  const auto& operator*() const { return *map_.find(*it_); }

  This &operator++() {
    ++it_;
    return *this;
  }

  bool operator==(const This &that) const { return it_ == that.it_; }

 private:
  const Map& map_;
  Iterator it_;
};

template <typename Map, typename Iterator>
auto make_key_to_key_value_iterator(const Map &map, const Iterator &it) {
  return Key2KeyValueIterator<Map, Iterator>(map, it);
}

template <typename Map>
class StateMap : public ScoresFilterMapMixin<Map> {
 public:
  auto states() const {
    const Map *p = crtp_this();
    return llvm::make_range(make_state_iterator(*p, p->cbegin()),
                            make_state_iterator(*p, p->cend()));
  }

  template <typename Container>
  auto states(const Container &container) const {
    const Map *p = crtp_this();
    auto b = make_key_to_key_value_iterator(*p, container.cbegin());
    auto e = make_key_to_key_value_iterator(*p, container.cend());
    return llvm::make_range(make_state_iterator(*p, b),
                            make_state_iterator(*p, e));
  }
  // auto states() const {  // TODO implement this as a generator
  //   using GraphCursor = typename Map::key_type;
  //   std::vector<State<GraphCursor>> states;
  //   for (const auto &kv : *crtp_this()) {
  //     states.push_back(Map::kv2state(kv));
  //   }
  //
  //   return states;
  // }

  // template <typename Container>
  // auto states(const Container &container) const {  // TODO implement this as a generator
  //   using GraphCursor = typename Map::key_type;
  //   std::vector<State<GraphCursor>> states;
  //   for (const auto &cursor : container) {
  //     auto it = crtp_this()->find(cursor);
  //     assert(it != crtp_this()->cend());
  //     states.push_back(Map::kv2state(*it));
  //   }
  //
  //   return states;
  // }

 private:
  Map *crtp_this() { return static_cast<Map *>(this); }
  const Map *crtp_this() const { return static_cast<const Map *>(this); }
};

template <typename GraphCursor>
class DeletionStateSet : public std::unordered_map<GraphCursor, ScoredPLink<GraphCursor>>,
                         public StateMap<DeletionStateSet<GraphCursor>> {
 public:
  template <typename KV>
  static State<GraphCursor> kv2state(const KV &kv) {
    const auto &cursor = kv.first;
    const auto &score = kv.second.score;
    const auto &plink = kv.second.plink;

    return {cursor, plink, score};
  }

  bool update(const GraphCursor &cursor,
              score_t score,
              const PathLinkRef<GraphCursor> &plink) {
    auto it_fl = this->insert({cursor, {plink, score}});
    auto it = it_fl.first;
    bool inserted = it_fl.second;

    if (inserted) {
      return true;
    }

    const score_t &prev_score = it->second.score;
    if (prev_score > score) {
      it->second = {plink, score};
      return true;
    }
    return false;
  }

  template <typename Set>
  size_t merge(const Set &S, score_t fee = 0) {
    size_t count = 0;

    for (const auto &state : S.states()) {
      count += update(state.cursor, state.score + fee, state.plink);
    }

    return count;
  }

  void increment(score_t fee = 0) {
    for (auto &kv : *this) {
      kv.second.score += fee;
    }
  }

  template <typename KV>  // TODO Use exact type here instof duck typing
  static score_t score_fnc(const KV &kv) {
    return kv.second.score;
  }
};

template <typename GraphCursor>
class StateSet : public std::unordered_map<GraphCursor, PathLinkRef<GraphCursor>>,
                 public StateMap<StateSet<GraphCursor>> {
 public:
  template <typename KV>
  static State<GraphCursor> kv2state(const KV &kv) {
    const auto &cursor = kv.first;
    const auto &score = kv.second->score();
    const auto &plink = kv.second;

    return {cursor, plink, score};
  }

  void set_event(size_t m, EventType type) {
    for (auto &kv : *this) {
      if (!kv.first.is_empty()) {
        kv.second->set_emission(m, type);
      }
    }
  }

  StateSet clone() const {
    StateSet copy{*this};
    for (auto &kv : copy) {
      kv.second = kv.second->clone();
    }

    return copy;
  }

  bool update(const GraphCursor &cursor, score_t score, const GraphCursor &from,
              const PathLinkRef<GraphCursor> &plink) {
    auto it_fl = this->insert({cursor, nullptr});
    auto it = it_fl.first;
    bool inserted = it_fl.second;
    if (inserted) {
      it->second = PathLink<GraphCursor>::create();
    }

    score_t prev = inserted ? std::numeric_limits<score_t>::infinity() : it->second->score();
    it->second->update(from, score, plink);
    return prev > score;
  }

  template <typename KV>  // TODO Use exact type here instof duck typing
  static score_t score_fnc(const KV &kv) {
    return kv.second->score();
  }

};

template <typename GraphCursor>
PathSet<GraphCursor> find_best_path(const hmm::Fees &fees,
                                    const std::vector<GraphCursor> &initial_original,
                                    typename GraphCursor::Context context) {
  using StateSet = StateSet<GraphCursor>;
  using DeletionStateSet = DeletionStateSet<GraphCursor>;
  const auto &code = fees.code;

  INFO("pHMM size: " << fees.M);
  if (!fees.check_i_loop(0)) {
    WARN("Positive-score insertion at the beginning");
  }
  if (!fees.check_i_loop(fees.M)) {
    WARN("Positive-score insertion at the end");
  }

  for (size_t i = 0; i <= fees.M; ++i) {
    if (!fees.check_i_loop(i)) {
      WARN("Positive-score insertion at position " << i);
    }
  }

  if (!fees.check_i_negative_loops()) {
    WARN("MODEL CONTAINS POSITIVE-SCORE I-LOOPS");
  }

  std::vector<GraphCursor> initial;

  auto transfer = [&code, &initial,context](StateSet &to, const auto &from, double transfer_fee,
                                            const std::vector<double> &emission_fees) {
    DEBUG_ASSERT((void*)(&to) != (void*)(&from), hmmpath_assert{});
    for (const auto &state : from.states()) {
      for (const auto &next : (state.cursor.is_empty() ? initial : state.cursor.next(context))) {
        double cost = state.score + transfer_fee + emission_fees[code(next.letter(context))];
        to.update(next, cost, state.cursor, state.plink);
      }
    }
  };

  auto transfer_upd = [&code,context](StateSet &to, const StateSet &from, double transfer_fee,
                                      const std::vector<double> &emission_fees, const std::string &,
                                      const std::unordered_set<GraphCursor> &keys) {
    DEBUG_ASSERT(&to != &from, hmmpath_assert{});
    std::unordered_set<GraphCursor> updated;
    for (const auto &state : from.states(keys)) {
      const auto &cur = state.cursor;
      const auto &fee = state.score;
      const auto &id = state.plink;
      for (const auto &next : cur.next(context)) {
        char letter = next.letter(context);
        double cost = fee + transfer_fee + emission_fees[code(letter)];
        if (to.update(next, cost, cur, id)) {
          updated.insert(next);
        }
      }
    }
    return updated;
  };

  auto i_loop_processing_negative = [&transfer_upd, &fees](StateSet &I, size_t m) {
    // TODO Review it
    const size_t max_insertions = 30;

    std::unordered_set<GraphCursor> updated;
    for (const auto &kv : I) {
      updated.insert(kv.first);
    }
    I.set_event(m, EventType::INSERTION);
    StateSet Inew = I.clone();
    for (size_t i = 0; i < max_insertions; ++i) {
      updated = transfer_upd(Inew, I, fees.t[m][p7H_II], fees.ins[m], "o", updated);
      Inew.set_event(m, EventType::INSERTION);
      TRACE(updated.size() << " items updated");
      for (const auto &cur : updated) {
        I[cur] = Inew[cur]->clone();  // TODO Implement minor updation detection
      }
    }
    I = std::move(Inew);  // It is necessary to copy minorly updated states
  };

  size_t not_cool_global_n_const = size_t(-1);
  auto i_loop_processing_non_negative_formally_correct_but_slow_and_potentially_leaking = [&fees, &code, context, &not_cool_global_n_const](StateSet &I, size_t m, const auto &filter) {
    const auto &emission_fees = fees.ins[m];
    const auto &transfer_fee = fees.t[m][p7H_II];

    TRACE(I.size() << " I states initially present in I-loop m = " << m);

    struct QueueElement {
      GraphCursor current_cursor;
      double score;

      bool operator<(const QueueElement &other) const {
        return this->score > other.score;
      }
    };

    std::priority_queue<QueueElement> q;

    for (const auto &kv : I) {
      const auto &current_cursor = kv.first;
      const score_t &score = kv.second->score();
      if (score > fees.absolute_threshold) {
        continue;
      }
      if (!filter(current_cursor)) {
        q.push({current_cursor, score});
      }
    }
    TRACE(q.size() << " I values in queue m = " << m);

    std::unordered_set<GraphCursor> processed;
    size_t taken_values = 0;
    while (!q.empty() && processed.size() < not_cool_global_n_const) {
      QueueElement elt = q.top();
      q.pop();
      ++taken_values;

      if (elt.score > fees.absolute_threshold) {
        break;
      }
      auto it_fl = processed.insert(elt.current_cursor);
      if (!it_fl.second) { // Already there
        continue;
      }

      const auto &id = I[elt.current_cursor];
      for (const auto &next : elt.current_cursor.next(context)) {
        char letter = next.letter(context);
        double cost = elt.score + transfer_fee + emission_fees[code(letter)];
        if (!filter(next)) {
          bool updated = I.update(next, cost, elt.current_cursor, id);
          // FIXME potential memory leak here!
          if (updated) {
            q.push({next, cost});
          }
        }
      }
    }

    TRACE(processed.size() << " states processed in I-loop m = " << m);
    TRACE(taken_values << " values extracted from queue m = " << m);
    // TODO update secondary references.
    // Cycle references may appear =(
  };

  auto i_loop_processing_non_negative = [&fees, &code, context](StateSet &I, size_t m, const auto &filter) {
    const auto &emission_fees = fees.ins[m];
    const auto &transfer_fee = fees.t[m][p7H_II];

    TRACE(I.size() << " I states initially present in I-loop m = " << m);
    std::unordered_set<GraphCursor> updated;

    struct QueueElement {
      GraphCursor current_cursor;
      double score;
      GraphCursor source_cursor;
      PathLinkRef<GraphCursor> source_state;

      bool operator<(const QueueElement &other) const {
        return this->score > other.score;
      }
    };

    std::priority_queue<QueueElement> q;

    for (const auto &kv : I) {
      const auto &current_cursor = kv.first;
      auto best = kv.second->best_ancestor();
      const auto &score = best->second.first;
      if (score > fees.absolute_threshold) {
        continue;
      }
      if (!filter(current_cursor)) {
        q.push({current_cursor, score, best->first, best->second.second});
      }
    }
    TRACE(q.size() << " I values in queue m = " << m);

    std::unordered_set<GraphCursor> processed;
    size_t taken_values = 0;
    while(!q.empty()) {
      QueueElement elt = q.top();
      q.pop();
      ++taken_values;

      if (elt.score > fees.absolute_threshold) {
        break;
      }

      auto it_fl = processed.insert(elt.current_cursor);
      if (!it_fl.second) { // Already there
        continue;
      }

      I.update(elt.current_cursor, elt.score, elt.source_cursor, elt.source_state);  // TODO return iterator to inserted/updated elt
      const auto &id = I[elt.current_cursor];
      for (const auto &next : elt.current_cursor.next(context)) {
        if (processed.count(next)) {
          continue;
        }
        char letter = next.letter(context);
        double cost = elt.score + transfer_fee + emission_fees[code(letter)];
        if (!filter(next)) {
          q.push({next, cost, elt.current_cursor, id});
        }
      }
    }

    TRACE(processed.size() << " states processed in I-loop m = " << m);
    TRACE(taken_values << " values extracted from queue m = " << m);
    // TODO update secondary references.
    // Cycle references may appear =(
  };

  auto i_loop_processing = [&](StateSet &I, size_t m, const auto &filter) {
    if (fees.is_i_loop_non_negative(m)) {
      if (fees.use_experimental_i_loop_processing) {
        return i_loop_processing_non_negative_formally_correct_but_slow_and_potentially_leaking(I, m, filter);
      } else {
        return i_loop_processing_non_negative(I, m, filter);
      }
    } else {
      return i_loop_processing_negative(I, m);
    }
  };

  auto dm_new = [&](DeletionStateSet &D, StateSet &M, const StateSet &I, size_t m) {
    DeletionStateSet preM = D;

    D.increment(fees.t[m - 1][p7H_DD]);
    D.merge(M, fees.t[m - 1][p7H_MD]);

    preM.increment(fees.t[m - 1][p7H_DM]);
    preM.merge(M, fees.t[m - 1][p7H_MM]);
    preM.merge(I, fees.t[m - 1][p7H_IM]);

    M.clear();
    transfer(M, preM, 0, fees.mat[m]);
  };

  depth_filter::impl::DepthAtLeast<GraphCursor> depth;

  INFO("Original (before filtering) initial set size: " << initial_original.size());
  // depth_filter::impl::Depth<GraphCursor> depth_naive;
  // for (const auto &cursor : initial_original) {
  //   INFO("Cursor " << cursor << "Depth: " << depth_naive.depth(cursor, context));
  // }
  double required_cursor_depth = static_cast<double>(fees.M) * fees.depth_filter_rate - fees.depth_filter_constant;
  if (fees.local) {
    required_cursor_depth = fees.minimal_match_length;
  }
  INFO("Depth required: " << required_cursor_depth);
  std::copy_if(initial_original.cbegin(), initial_original.cend(), std::back_inserter(initial),
               [&](const GraphCursor &cursor) { return depth.depth_at_least(cursor, required_cursor_depth,
                                                                            context); });
  INFO("Initial set size: " << initial.size());

  StateSet I, M;
  DeletionStateSet D;
  auto source = PathLink<GraphCursor>::create_source();
  auto sink = PathLink<GraphCursor>::create();
  M[GraphCursor()] = source;
  auto update_sink = [&sink](const auto &S, double fee) {
    for (const auto &state : S.states()) {
      sink->update(state.cursor, state.score + fee, state.plink);
    }
  };

  INFO("The number of links (M): " << fees.M);

  size_t positions_left = fees.M;
  auto depth_filter_cursor = [&](const GraphCursor &cursor) -> bool {
    double required_cursor_depth = static_cast<double>(positions_left) * fees.depth_filter_rate - fees.depth_filter_constant;
    if (fees.local) {
      required_cursor_depth = 0;
    }
    return !depth.depth_at_least(cursor, required_cursor_depth, context);
  };

  transfer(I, M, fees.t[0][p7H_MI], fees.ins[0]);
  i_loop_processing(I, 0, depth_filter_cursor);  // Do we really need I at the beginning???
  I.set_event(0, EventType::INSERTION);
  size_t n = 1;
  for (size_t m = 1; m <= fees.M; ++m) {
    positions_left = fees.M - m;

    if (fees.local && m > 1) {
      D.update(GraphCursor(), fees.cleavage_cost, source);
    }
    dm_new(D, M, I, m);

    I.clear();
    transfer(I, M, fees.t[m][p7H_MI], fees.ins[m]);
    i_loop_processing(I, m, depth_filter_cursor);

    size_t n_of_states = D.size() + I.size() + M.size();

    TRACE("# states " << m << " => " << n_of_states);
    size_t top = n_of_states;
    if (m > 25) {
      not_cool_global_n_const = top = fees.state_limits.l25;
    }
    if (m > 100) {
      not_cool_global_n_const = top = fees.state_limits.l100;
    }
    if (m > 500) {
      not_cool_global_n_const = top = fees.state_limits.l500;
    }

    if (m >= n) {
      INFO("Step #: " << m);
      INFO("# states " << m << " => " << n_of_states << ": I = " << I.size() << " M = " << M.size() << " D = " << D.size());
    }

    I.score_filter(top, fees.absolute_threshold);
    M.score_filter(top, fees.absolute_threshold);
    D.score_filter(top, fees.absolute_threshold);

    size_t depth_filtered = 0;
    if (m % 8 == 0) {
      depth_filtered += I.filter_key(depth_filter_cursor);
      depth_filtered += M.filter_key(depth_filter_cursor);
      depth_filtered += D.filter_key(depth_filter_cursor);
    }

    I.set_event(m, EventType::INSERTION);
    M.set_event(m, EventType::MATCH);

    // collapse and trim
    for (auto &kv : I) {
      kv.second->collapse_and_trim();
    }
    for (auto &kv : M) {
      kv.second->collapse_and_trim();
    }

    if (fees.local) {
      update_sink(D, fees.cleavage_cost);
      sink->collapse_and_trim();
    }

    if (m >= n) {
      INFO("depth-filtered " << depth_filtered << ", positions left = " << positions_left << " states m = " << m);
      INFO("I = " << I.size() << " M = " << M.size() << " D = " << D.size());
      auto scores = M.scores();
      std::sort(scores.begin(), scores.end());
      if (scores.size() > 100) { scores.resize(100); }
      for (auto &score : scores) {
        score = -score;
      }
      INFO("Top scores: " << scores);
      n <<= 1;
    }
  }

  INFO("Max stack size in Depth: " << depth.max_stack_size());

  PathSet<GraphCursor> result(sink);

  update_sink(D, fees.t[fees.M][p7H_DM]);
  update_sink(I, fees.t[fees.M][p7H_IM]);  // Do we really need I at the end?
  update_sink(M, fees.t[fees.M][p7H_MM]);
  sink->collapse_and_trim();
  // sink.apply([](auto &p) { const_cast<PathLink<GraphCursor>&>(p).collapse_and_trim(); });

  DEBUG(sink->object_count_current() << " pathlink objects");
  DEBUG(sink->object_count_max() << " pathlink objects maximum");
  DEBUG(sink->object_count_constructed() << " pathlink objects constructed");

  return result;
}

}  // namespace impl

// vim: set ts=2 sw=2 et :

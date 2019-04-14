#pragma once

#include "fees.hpp"
#include "pathtree.hpp"
#include "depth_filter.hpp"
#include "cursor_utils.hpp"

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

#include <parallel_hashmap/phmap.h>

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
class DeletionStateSet : public phmap::flat_hash_map<GraphCursor, ScoredPLink<GraphCursor>>,
                         public StateMap<DeletionStateSet<GraphCursor>> {
 public:
  size_t collapse_all() {
    size_t count = 0;
    for (auto &kv : *this) {
      count += kv.second.plink->collapse_and_trim();
    }
    return count;
  }

  template <typename KV>
  static State<GraphCursor> kv2state(const KV &kv) {
    const auto &cursor = kv.first;
    const auto &score = kv.second.score;
    const auto &plink = kv.second.plink;

    return {cursor, plink, score};
  }

  bool update(score_t score,
              const PathLinkRef<GraphCursor> &plink) {
    const GraphCursor &cursor = plink->cursor();
    auto it_fl = this->insert({cursor, {plink, score}});
    const auto &it = it_fl.first;
    const bool &inserted = it_fl.second;

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
      count += update(state.score + fee, state.plink);
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
class StateSet : public phmap::flat_hash_map<GraphCursor, PathLinkRef<GraphCursor>>,
                 public StateMap<StateSet<GraphCursor>> {
 public:
  size_t collapse_all() {
    size_t count = 0;
    for (auto &kv : *this) {
      count += kv.second->collapse_and_trim();
    }
    return count;
  }

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

  score_t get_cost(const GraphCursor &cursor) const {
    auto it = this->find(cursor);
    return it != this->cend() ? it->second->score() : std::numeric_limits<score_t>::infinity();
  }

  StateSet clone() const {
    StateSet copy{*this};
    for (auto &kv : copy) {
      kv.second = kv.second->clone();
    }

    return copy;
  }

  bool update(const GraphCursor &cursor, score_t score,
              const PathLinkRef<GraphCursor> &plink) {
    auto it_fl = this->insert({cursor, nullptr});
    const auto &it = it_fl.first;
    const bool &inserted = it_fl.second;
    if (inserted) {
      it->second = PathLink<GraphCursor>::create(cursor);
    }

    score_t prev = inserted ? std::numeric_limits<score_t>::infinity() : it->second->score();
    it->second->update(score, plink);
    return prev > score;
  }

  template <typename KV>  // TODO Use exact type here instof duck typing
  static score_t score_fnc(const KV &kv) {
    return kv.second->score();
  }

};

template <typename GraphCursor>
PathSet<GraphCursor> find_best_path(const hmm::Fees &fees,
                                    const std::vector<GraphCursor> &cursors,
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

  auto vcursors = vertex_cursors(cursors, context);
  INFO("Vertex cursors: " << vcursors.size() << "/" << cursors.size());

  depth_filter::DepthInt<GraphCursor> depth;

  std::vector<GraphCursor> initial;

  auto transfer = [&code, &initial,context](StateSet &to, const auto &from, double transfer_fee,
                                            const std::vector<double> &emission_fees) {
    DEBUG_ASSERT((void*)(&to) != (void*)(&from), hmmpath_assert{});
    for (const auto &state : from.states()) {
      for (const auto &next : (state.cursor.is_empty() ? initial : state.cursor.next(context))) {
        double cost = state.score + transfer_fee + emission_fees[code(next.letter(context))];
        to.update(next, cost, state.plink);
      }
    }
  };

  auto loop_transfer_ff= [&code, context, &fees, &depth, &vcursors](StateSet &I, double transfer_fee,
                                                                    const std::vector<double> &emission_fees,
                                                                    const phmap::flat_hash_set<GraphCursor> &keys) {
    DEBUG("loop_transfer_ff begins");
    StateSet Inext;
    std::vector<GraphCursor> updated_vertices;
    std::vector<GraphCursor> updated_nonvertices;
    phmap::flat_hash_set<GraphCursor> relaxed;

    std::vector<GraphCursor> stack = extract_leftmost_cursors(keys, vcursors, context);
    while (!stack.empty()) {
      GraphCursor cursor = stack.back(); stack.pop_back();
      VERIFY(I.count(cursor));
      const auto plink = I[cursor];  // Do not take reference here, we use flat map, reference would be invalidated after insertion!!!
      VERIFY(plink);
      relaxed.insert(cursor);

      double required_cursor_depth = static_cast<double>(fees.minimal_match_length) - 1 - static_cast<double>(plink->max_prefix_size());
      for (const auto &next : cursor.next(context)) {
        double cost = plink->score() + transfer_fee + emission_fees[code(next.letter(context))];
        if (!vcursors.count(next)) {
          VERIFY(next.prev(context).size() == 1 && next.next(context).size() == 1);
          DEBUG("FAST FORWARD");
          bool successful_update = (cost <= fees.absolute_threshold) &&
                                   (depth.depth_at_least(next, required_cursor_depth, context)) &&
                                   (cost < I.get_cost(next));
          if (successful_update) {
            I[next] = PathLink<GraphCursor>::create(next);
            I[next]->update(cost, plink);
            updated_nonvertices.push_back(next);
            stack.push_back(next);
            DEBUG("next " << next << " added into stack, stack size " << stack.size());
          } else {
            // Go to the next key
            DEBUG("Update inefficient; go forward along the edge to the next non-vertex key");
            GraphCursor nn = next;
            while (!vcursors.count(nn) && !keys.count(nn)) {
              nn = nn.next(context)[0];  // It's correct due to triviality of nn
            }
            if (keys.count(nn) && !vcursors.count(nn)) {
              stack.push_back(nn);
            }
          }
        } else {
          if (cost > fees.absolute_threshold) continue;
          if (!depth.depth_at_least(next, required_cursor_depth, context)) continue;
          Inext.update(next, cost, plink);
          if (cost < I.get_cost(next)) {
            updated_vertices.push_back(next);
          }
        }
      }
    }
    for (const GraphCursor &cursor : keys) {
      VERIFY(relaxed.count(cursor));
    }

    remove_duplicates(updated_vertices);
    // updated_nonvertices cannot contain duplicates TODO VERIFY it

    for (const GraphCursor &cursor : updated_vertices) {
      auto plink = std::move(Inext[cursor]);
      plink->collapse_and_trim();
      I[cursor] = std::move(plink);
    }

    phmap::flat_hash_set<GraphCursor> updated;
    updated.insert(updated_vertices.cbegin(), updated_vertices.cend());
    updated.insert(updated_nonvertices.cbegin(), updated_nonvertices.cend());

    return updated;
  };

  auto loop_transfer_negative = [&code, context, &fees, &depth](StateSet &I, double transfer_fee,
                                                                const std::vector<double> &emission_fees,
                                                                const auto &keys,
                                                                bool just_all = false) {
    StateSet Inext;
    std::vector<GraphCursor> updated;
    auto process = [&](const auto &collection) -> void {
      for (const auto &state : collection) {
        double required_cursor_depth = static_cast<double>(fees.minimal_match_length) - 1 - static_cast<double>(state.plink->max_prefix_size());
        for (const auto &next : state.cursor.next(context)) {
          double cost = state.score + transfer_fee + emission_fees[code(next.letter(context))];
          if (cost > fees.absolute_threshold) continue;
          if (!depth.depth_at_least(next, required_cursor_depth, context)) continue;
          Inext.update(next, cost, state.plink);
          if (cost < I.get_cost(next)) {
            updated.push_back(next);
          }
        }
      }
    };
    just_all ? process(I.states()) : process(I.states(keys));
    remove_duplicates(updated);

    for (const GraphCursor &cursor : updated) {
      auto plink = std::move(Inext[cursor]);
      plink->collapse_and_trim();
      I[cursor] = std::move(plink);
    }

    return updated;
  };

  auto i_loop_processing_ff = [&loop_transfer_ff, &fees](StateSet &I, size_t m) {
    phmap::flat_hash_set<GraphCursor> updated;
    for (const auto &kv : I) {
      updated.insert(kv.first);
    }
    I.set_event(m, EventType::INSERTION);
    for (size_t i = 0; i < fees.max_insertion_length && !updated.empty(); ++i) {
      updated = loop_transfer_ff(I, fees.t[m][p7H_II], fees.ins[m], updated);
      if (is_power_of_two_or_zero(m)) {
        INFO("Updated: " << updated.size() << " over " << I.size() << " on i = " << i << " m = " << m);
      }
      for (const GraphCursor &cursor : updated) {
        VERIFY(I.count(cursor));
        I[cursor]->set_emission(m, EventType::INSERTION);
      }
    }
  };

  auto i_loop_processing_universal = [&loop_transfer_negative, &fees](StateSet &I, size_t m) {
    std::vector<GraphCursor> updated;
    I.set_event(m, EventType::INSERTION);
    for (size_t i = 0; i < fees.max_insertion_length; ++i) {
      updated = loop_transfer_negative(I, fees.t[m][p7H_II], fees.ins[m], updated, /*just_all*/ i == 0);
      if (is_power_of_two_or_zero(m)) {
        INFO("Updated: " << updated.size() << " over " << I.size() << " on i = " << i << " m = " << m);
      }
      for (const GraphCursor &cursor : updated) {
        I[cursor]->set_emission(m, EventType::INSERTION);
      }
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

  INFO("Original (before filtering) initial set size: " << cursors.size());
  // depth_filter::impl::Depth<GraphCursor> depth_naive;
  // for (const auto &cursor : cursors) {
  //   INFO("Cursor " << cursor << "Depth: " << depth_naive.depth(cursor, context));
  // }
  double required_cursor_depth = static_cast<double>(fees.minimal_match_length);
  INFO("Depth required: " << required_cursor_depth);
  std::copy_if(cursors.cbegin(), cursors.cend(), std::back_inserter(initial),
               [&](const GraphCursor &cursor) { return depth.depth_at_least(cursor, required_cursor_depth,
                                                                            context); });
  INFO("Initial set size: " << initial.size());

  StateSet I, M;
  DeletionStateSet D;
  auto source = PathLink<GraphCursor>::create_source();
  auto sink = PathLink<GraphCursor>::create_sink();
  M[GraphCursor()] = source;
  auto update_sink = [&sink](const auto &S, double fee) {
    for (const auto &state : S.states()) {
      sink->update(state.score + fee, state.plink);
    }
  };

  INFO("The number of links (M): " << fees.M);

  auto depth_filter_kv = [&](const auto &cursor_value) -> bool {
    const GraphCursor &cursor = cursor_value.first;
    size_t prefix_len = cursor_value.second->max_prefix_size();
    double required_cursor_depth = static_cast<double>(fees.minimal_match_length) - static_cast<double>(prefix_len);
    // bool dal = depth.depth_at_least(cursor, required_cursor_depth, context);
    // bool di = depth_int.depth_at_least(cursor, required_cursor_depth, context);
    // if (di != dal) {
    //   INFO(di << " " << dal << " " << cursor << " " << required_cursor_depth << " " << depth_int.depth(cursor, context));
    // }
    // VERIFY(di == dal);
    return !depth.depth_at_least(cursor, required_cursor_depth, context);
  };

  transfer(I, M, fees.t[0][p7H_MI], fees.ins[0]);
  i_loop_processing_ff(I, 0);  // Do we really need I at the beginning???
  I.set_event(0, EventType::INSERTION);
  for (size_t m = 1; m <= fees.M; ++m) {
    if (fees.local && m > 1) {  // FIXME check latter condition. Does it really make sense?
      D.update(fees.cleavage_cost, source);
    }
    dm_new(D, M, I, m);
    M.collapse_all();

    I.clear();
    transfer(I, M, fees.t[m][p7H_MI], fees.ins[m]);
    I.collapse_all();
    i_loop_processing_ff(I, m);

    size_t n_of_states = D.size() + I.size() + M.size();

    TRACE("# states " << m << " => " << n_of_states);
    size_t top = n_of_states;
    if (m > 25) {
      top = fees.state_limits.l25;
    }
    if (m > 100) {
      top = fees.state_limits.l100;
    }
    if (m > 500) {
      top = fees.state_limits.l500;
    }

    if (is_power_of_two_or_zero(m)) {
      INFO("Step #: " << m);
      INFO("# states " << m << " => " << n_of_states << ": I = " << I.size() << " M = " << M.size() << " D = " << D.size());
    }

    I.score_filter(top, fees.absolute_threshold);
    M.score_filter(top, fees.absolute_threshold);
    D.score_filter(top, fees.absolute_threshold);

    size_t depth_filtered = 0;
    if (m % 1 == 0) {
      depth_filtered += I.filter_key_value(depth_filter_kv);
      depth_filtered += M.filter_key_value(depth_filter_kv);
      // depth_filtered += D.filter_key_value(depth_filter_kv);  // depth filter for Ds is not required
    }

    I.set_event(m, EventType::INSERTION);
    M.set_event(m, EventType::MATCH);

    if (fees.local) {
      update_sink(D, fees.cleavage_cost);
      sink->collapse_and_trim();  // FIXME subtract cost for transition -> D state ?
    }

    if (is_power_of_two_or_zero(m)) {
      INFO("depth-filtered " << depth_filtered << " position in HMM " << m);
      INFO("I = " << I.size() << " M = " << M.size() << " D = " << D.size());
      auto scores = M.scores();
      std::sort(scores.begin(), scores.end());
      if (scores.size() > 100) { scores.resize(100); }
      for (auto &score : scores) {
        score = -score;
      }
      INFO("Top scores: " << scores);
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

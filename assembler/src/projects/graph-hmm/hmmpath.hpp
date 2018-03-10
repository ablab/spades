#pragma once

#include "cursor.hpp"
#include "fees.hpp"
#include "pathtree.hpp"

#include <unordered_map>
#include <memory>
#include <vector>

namespace impl {

using pathtree::PathLink;

template <typename GraphCursor>
class StateSet : public std::unordered_map<GraphCursor, std::shared_ptr<PathLink<GraphCursor>>> {
 public:
  const std::shared_ptr<PathLink<GraphCursor>> &get_or_create(const GraphCursor &key) {
    return this->insert({key, std::make_shared<PathLink<GraphCursor>>()}).first->second;
  }

  StateSet clone() const {
    StateSet copy{*this};
    for (auto &kv : copy) {
      kv.second = kv.second->clone();
    }

    return copy;
  }

  bool update(const GraphCursor &key, double score, GraphCursor from,
              const std::shared_ptr<PathLink<GraphCursor>> &traj) {
    auto it_fl = this->insert({key, std::make_shared<PathLink<GraphCursor>>()});
    double prev = it_fl.second ? std::numeric_limits<double>::infinity() : it_fl.first->second->score();
    // double prev = it_fl.first->second->score();
    it_fl.first->second->update(from, score, traj);
    return prev > score;
  }

  // TODO implement method detecting upd of any kind
};

template <typename GraphCursor>
StateSet<GraphCursor> top_filter(const StateSet<GraphCursor> &S, size_t top, double threshold) {
  using StateSet = StateSet<GraphCursor>;
  top = std::min(top, S.size());
  std::vector<std::pair<typename StateSet::key_type, typename StateSet::mapped_type>> v(S.cbegin(), S.cend());

  // TODO Use std::nth_element
  std::sort(v.begin(), v.end(), [](const auto &e1, const auto &e2) { return e1.second->score() < e2.second->score(); });
  if (v.size() < top) {
    top = v.size();
  }
  size_t size = 0;
  while (size < top && v[size].second->score() < threshold) {
    ++size;
  }
  v.resize(size);

  StateSet result;
  result.insert(CONST_ALL(v));
  return result;
}

template <typename GraphCursor>
std::vector<std::pair<std::string, double>> find_best_path(const hmm::Fees &fees, const std::vector<GraphCursor> &initial) {
  using StateSet = StateSet<GraphCursor>;
  const auto &code = fees.code;

  INFO("pHMM size: " << fees.M);
  if (!fees.check_i_loop(0)) {
    WARN("Negative-cost insertion at the beginning");
  }
  if (!fees.check_i_loop(fees.M)) {
    WARN("Negative-cost insertion at the end");
  }

  for (size_t i = 0; i <= fees.M; ++i) {
    if (!fees.check_i_loop(i)) {
      WARN("Negative-cost insertion at position " << i);
    }
  }

  if (!fees.check_i_negative_loops()) {
    WARN("MODEL CONTAINS NEGATIVE I-LOOPS");
  }

  auto transfer = [&code, &initial](StateSet &to, const StateSet &from, double transfer_fee,
                                    const std::vector<double> &emission_fees, const std::string & = "") {
    assert(&to != &from);
    for (const auto &kv : from) {
      const auto &cur = kv.first;
      const auto &fee = kv.second->score();
      const auto &id = kv.second;
      if (cur.is_empty()) {
        // This branch is used only during BEGIN->M, BEGIN->I and D->M transfers
        for (size_t i = 0; i < initial.size(); ++i) {
          const auto &next = initial[i];
          double cost = fee + transfer_fee + emission_fees[code(next.letter())];
          to.update(next, cost, cur, id);
        }
      } else {
        auto next_pairs = cur.next_pairs();
        for (size_t i = 0; i < next_pairs.size(); ++i) {
          const auto &next = next_pairs[i].first;
          char letter = next_pairs[i].second;
          double cost = fee + transfer_fee + emission_fees[code(letter)];
          to.update(next, cost, cur, id);
        }
      }
    }
  };
  auto transfer_upd = [&code, &initial](StateSet &to, const StateSet &from, double transfer_fee,
                                        const std::vector<double> &emission_fees, const std::string &,
                                        const std::unordered_set<GraphCursor> &keys) {
    assert(&to != &from);
    std::unordered_set<GraphCursor> updated;
    for (const auto &cur : keys) {
      auto it = from.find(cur);
      assert(it != from.cend());
      const auto &fee = it->second->score();
      const auto &id = it->second;
      auto next_pairs = cur.next_pairs();
      for (size_t i = 0; i < next_pairs.size(); ++i) {
        const auto &next = next_pairs[i].first;
        char letter = next_pairs[i].second;
        double cost = fee + transfer_fee + emission_fees[code(letter)];
        if (to.update(next, cost, cur, id)) {
          updated.insert(next);
        }
      }
    }
    return updated;
  };

  auto i_loop_processing = [&transfer_upd, &fees](StateSet &I, size_t m) {
    const size_t max_insertions = 30;

    std::unordered_set<GraphCursor> updated;
    for (const auto &kv : I) {
      updated.insert(kv.first);
    }
    StateSet Inew = I.clone();
    for (size_t i = 0; i < max_insertions; ++i) {
      updated = transfer_upd(Inew, I, fees.t[m][p7H_II], fees.ins[m], "o", updated);
      TRACE(updated.size() << " items updated");
      for (const auto &cur : updated) {
        I[cur] = Inew[cur]->clone();  // TODO Implement minor updation detection
      }
    }
    I = std::move(Inew);  // It is necessary to copy minorly updated states
  };

  auto merge_state_set = [](StateSet &target, const StateSet &source, double transfer_fee = 0) {
    for (const auto &kv : source) {
      const auto &cur = kv.first;
      const auto &id = kv.second;
      target.get_or_create(cur)->merge_update(id.get(), transfer_fee);
    }
  };

  auto dm_new = [&](StateSet &D, StateSet &M, const StateSet &I, size_t m) {
    StateSet Dnew;
    merge_state_set(Dnew, M, fees.t[m - 1][p7H_MD]);
    merge_state_set(Dnew, D, fees.t[m - 1][p7H_DD]);

    StateSet Mnew;
    transfer(Mnew, M, fees.t[m - 1][p7H_MM], fees.mat[m], "m");
    transfer(Mnew, D, fees.t[m - 1][p7H_DM], fees.mat[m], "m");
    transfer(Mnew, I, fees.t[m - 1][p7H_IM], fees.mat[m], "m");

    M = std::move(Mnew);
    D = std::move(Dnew);
  };

  INFO("Initial set size: " << initial.size());

  StateSet I, M, D;
  const auto empty = GraphCursor();
  auto base = PathLink<GraphCursor>::master_source();
  M[empty] = base;  // TODO Implement and use empty Trajectory() instead of Trajectory(0)

  INFO("The number of links (M): " << fees.M);

  transfer(I, M, fees.t[0][p7H_MI], fees.ins[0], "i");
  i_loop_processing(I, 0);  // Do we really need I at the beginning???
  for (size_t m = 1; m <= fees.M; ++m) {
    INFO("Step #: " << m);

    dm_new(D, M, I, m);
    I.clear();
    transfer(I, M, fees.t[m][p7H_MI], fees.ins[m], "i");
    i_loop_processing(I, m);

    size_t top = std::max({D.size(), I.size(), M.size()});
    if (m > 10) {
      top = std::min<size_t>(1000000, top);
    }
    if (m > 50) {
      top = std::min<size_t>(10000, top);
    }
    if (m > 500) {
      top = std::min<size_t>(10000, top);
    }

    I = top_filter(I, top, 10);
    M = top_filter(M, top, 10);
    D = top_filter(D, top, 10);
  }

  PathLink<GraphCursor> terminal;
  auto upd_terminal = [&](const StateSet &S, double fee) {
    for (const auto &kv : S) {
      terminal.update(kv.first, kv.second->score() + fee, kv.second);
    }
  };

  upd_terminal(D, fees.t[fees.M][p7H_DM]);
  upd_terminal(I, fees.t[fees.M][p7H_DM]);  // Do we really need I at the end?
  upd_terminal(M, fees.t[fees.M][p7H_MM]);

  INFO("Best score: " << terminal.score());
  INFO("Best of the best");
  INFO(terminal.best_path_string());

  auto result = terminal.top_k_string(10000);
  return result;
}

}  // namespace impl

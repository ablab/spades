  // depth_filter::impl::DepthAtLeast<GraphCursor> depth;
  //
  // if (true) {
  //   depth_filter::impl::Depth<GraphCursor> depth_recursive;
  //   for (const auto &cursor : initial_original) {
  //     size_t di = depth.depth(cursor, context);
  //     double dd = depth_recursive.depth(cursor, context);
  //     size_t dd2di = dd == std::numeric_limits<double>::infinity() ? std::numeric_limits<size_t>::max() : static_cast<size_t>(dd);
  //     if (di != dd2di) {
  //       INFO(di << " " << dd2di);
  //     }
  //     if (!depth_at_least.depth_at_least(cursor, di, context)) {
  //       INFO("! depth_at_least");
  //     }
  //
  //     if (di != depth.INF && depth_at_least.depth_at_least(cursor, di + 1, context)) {
  //       INFO("depth_at_least + 1");
  //     }
  //   }
  // }

  auto transfer_upd = [&code,context](StateSet &to, const StateSet &from, double transfer_fee,
                                      const std::vector<double> &emission_fees,
                                      const std::unordered_set<GraphCursor> &keys) {
    DEBUG_ASSERT(&to != &from, hmmpath_assert{});
    std::unordered_set<GraphCursor> updated;
    for (const auto &state : from.states(keys)) {
      for (const auto &next : state.cursor.next(context)) {
        double cost = state.score + transfer_fee + emission_fees[code(next.letter(context))];
        if (to.update(next, cost, state.plink)) {
          updated.insert(next);
        }
      }
    }
    return updated;
  };


  auto i_loop_processing_negative_old = [&transfer_upd, &fees](StateSet &I, size_t m) {
    // TODO Review it
    // FIXME add proper collapsing
    const size_t max_insertions = 30;

    std::unordered_set<GraphCursor> updated;
    for (const auto &kv : I) {
      updated.insert(kv.first);
    }
    I.set_event(m, EventType::INSERTION);
    StateSet Inew = I.clone();
    for (size_t i = 0; i < max_insertions; ++i) {
      updated = transfer_upd(Inew, I, fees.t[m][p7H_II], fees.ins[m], updated);
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
    DEBUG("Experimental I-loops processing called");
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
      if (!filter(kv)) {
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
        // if (!filter(next)) {
          bool updated = I.update(next, cost, id);
          // FIXME potential memory leak here!
          if (updated) {
            q.push({next, cost});
          }
        // }
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
      const auto &score = best->first;
      if (score > fees.absolute_threshold) {
        continue;
      }
      if (!filter(kv)) {
        q.push({current_cursor, score, best->second->cursor(), best->second});
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

      I.update(elt.current_cursor, elt.score, elt.source_state);  // TODO return iterator to inserted/updated elt
      const auto &id = I[elt.current_cursor];
      for (const auto &next : elt.current_cursor.next(context)) {
        if (processed.count(next)) {
          continue;
        }
        char letter = next.letter(context);
        double cost = elt.score + transfer_fee + emission_fees[code(letter)];
        // if (!filter(next)) {
          q.push({next, cost, elt.current_cursor, id});
        // }
      }
    }

    TRACE(processed.size() << " states processed in I-loop m = " << m);
    TRACE(taken_values << " values extracted from queue m = " << m);
    // TODO update secondary references.
    // Cycle references may appear =(
  };

  auto i_loop_processing = [&](StateSet &I, size_t m, const auto &filter) {
    return i_loop_processing_negative2(I, m);

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


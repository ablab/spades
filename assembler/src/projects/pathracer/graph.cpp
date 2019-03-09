#include "graph.hpp"
#include "hmmpath.hpp"

std::vector<DBGraph::GraphCursor> DBGraph::all() const {
  std::unordered_set<GraphCursor> result;

  for (size_t id = 0; id < edges_.size(); ++id) {
    const auto &edge = edges_[id];
    for (size_t i = 0; i < edge.size(); ++i) {
      result.insert(get_pointer(id, i));
    }
  }

  return std::vector<GraphCursor>(result.cbegin(), result.cend());
}

PathSet<ReversedGraphCursor<Graph::GraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                    const std::vector<ReversedGraphCursor<Graph::GraphCursor>> &initial,
                                                                    const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<ReversedGraphCursor<DBGraph::GraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                      const std::vector<ReversedGraphCursor<DBGraph::GraphCursor>> &initial,
                                                                      const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<DBGraph::GraphCursor> find_best_path(const hmm::Fees &fees,
                                             const std::vector<DBGraph::GraphCursor> &initial,
                                             const void *context) {
  return impl::find_best_path(fees, initial, context);
}

PathSet<Graph::GraphCursor> find_best_path(const hmm::Fees &fees,
                                           const std::vector<Graph::GraphCursor> &initial,
                                           const void *context) {
  return impl::find_best_path(fees, initial, context);
}

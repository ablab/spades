//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include <vector>

#include "coordinates_handler.hpp"
//#include "coloring.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"

namespace cap {

/*
 * ConditionedSmartSetIterator acts much like SmartSetIterator, but, unlike the above case,
 * one can (and must) provide merge handler that will decide whether to add merged edge to
 * the set being iterated or not (extending add_new_ parameter logic of SmartIterator)
 * Also has the ability to be `reset` (i.e. start from the begin-iterator with respect to
 * added and deleted values)
 * MergeHandler class/struct must provide:
 *  bool operator()(const std::vector<ElementId> &, ElementId)
 */
template<class Graph, typename ElementId, class MergeHandler>
class ConditionedSmartSetIterator : public SmartSetIterator<Graph, ElementId> {
    typedef SmartSetIterator<Graph, ElementId> base;

    MergeHandler &merge_handler_;
    std::unordered_set<ElementId> true_elements_;

  public:

    template <class Iterator>
    ConditionedSmartSetIterator(const Graph &graph, Iterator begin, Iterator end,
                                MergeHandler &merge_handler)
            : SmartSetIterator<Graph, ElementId>(graph, begin, end),
              merge_handler_(merge_handler),
              true_elements_() {

        for (auto it = begin; it != end; ++it) {
            true_elements_.insert(*it);
        }
    }

    void HandleAdd(ElementId v) override {
        TRACE("handleAdd " << this->g().str(v));
        if (true_elements_.count(v)) {
            this->push(v);
        }
    }

    void HandleDelete(ElementId v) override {
        TRACE("handleDel " << this->g().str(v));
        base::HandleDelete(v);
        true_elements_.erase(v);
    }

    void HandleMerge(const std::vector<ElementId>& old_edges, ElementId new_edge) override {
        TRACE("handleMer " << this->g().str(new_edge));
        if (merge_handler_(old_edges, new_edge)) {
            true_elements_.insert(new_edge);
        }
    }

private:
    DECL_LOGGER("ConditionedSmartSetIterator");
};

template <class Graph>
class PathProjector {
 public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
  typedef std::vector<EdgeId> Path;
  typedef typename CoordinatesHandler<Graph>::PosArray PosArray;
  typedef unsigned uint;

  PathProjector(Graph &g, CoordinatesHandler<Graph> &coordinates_handler)
      : g_(g),
        coordinates_handler_(coordinates_handler),
        edge_remover_(g),
        is_deleting_locked_(false) {
  }

  virtual std::vector<Path> FilterPaths(const std::vector<Path> &paths) const {
    return paths;
  }
  virtual std::vector<PosArray> GetThreadsToDelete(
      const std::vector<Path> &paths) const {
    std::vector<PosArray> threads_to_delete;
    for (const auto &path : paths) {
      threads_to_delete.push_back(coordinates_handler_.GetContiguousThreads(path));
    }
    return threads_to_delete;
  }

  virtual size_t ChooseBasePath(const std::vector<Path> &paths,
      const std::vector<PosArray> &threads_to_delete) const {
    std::vector<size_t> num_bridges;
    std::vector<size_t> sum_multiplicities;
    for (const auto &path : paths) {
      num_bridges.push_back(CalcBridges(path, threads_to_delete.back().size()));
      sum_multiplicities.push_back(CalcMultiplicitySum(path));
    }
    
    size_t chosen_path = 0;
    for (size_t i = 1; i < paths.size(); ++i) {
      if (num_bridges[i] < num_bridges[chosen_path] ||
            (num_bridges[i] == num_bridges[chosen_path] &&
            sum_multiplicities[i] > sum_multiplicities[chosen_path])) {
        chosen_path = i;
      }
    }

    return chosen_path;
  }

  bool CollapsePaths(const std::vector<Path> &paths_to_collapse) {
    TRACE("CollapsePaths Begin");

    const std::vector<Path> &paths = FilterPaths(paths_to_collapse);
    if (CheckPresenceOfSelfContiguousEdges(paths))
      return false;

    std::vector<PosArray> threads_to_delete =
        GetThreadsToDelete(paths);
    size_t chosen_path = ChooseBasePath(paths, threads_to_delete);

    // RC paths
    std::vector<Path> rc_paths;
    for (const auto &path : paths) {
      Path rc_path;
      for (auto it = path.rbegin(); it != path.rend(); ++it) {
        rc_path.push_back(g_.conjugate(*it));
      }
      rc_paths.push_back(rc_path);
    }
    auto rc_threads = GetThreadsToDelete(rc_paths);

    // We restrict merging RC and original genome threads
    if (!CheckForRCAbsence(threads_to_delete))
      return false;
    if (!CheckForCorrectPaths(paths))
      return false;


    if (!CheckDeletionOfIntouchables(paths, threads_to_delete, chosen_path))
      return false;
    std::unordered_set<size_t> bad_paths = GetIntersectingPaths(paths,
        chosen_path);

    if(false && bad_paths.size() == paths.size() - 1) {
    for (size_t i = 0; i < paths.size(); ++i) {
      INFO("failed: paths " << Debug(paths[i]) << " and " << Debug(rc_paths[i]));
      for (const auto &e : threads_to_delete[i])
        INFO("to delete 1: " << int(e.first) << " - " << e.second);
      for (const auto &e : rc_threads[i])
        INFO("to delete 2: " << int(e.first) << " - " << e.second);
      INFO("---------paths1:");
      for (const auto &e : paths[i])
        coordinates_handler_.DebugOutput(e);
      INFO("---------paths2:");
      for (const auto &e : rc_paths[i])
        coordinates_handler_.DebugOutput(e);
    }
      VERIFY(false);
    }

    DEBUG("Collapsing paths:");
    for (size_t i = 0; i < paths.size(); ++i) {
        if (i == chosen_path)
            continue;

        DEBUG(debug::Debug(g_, paths[i]));
        for (const auto &t : threads_to_delete[i]) {
            DEBUG("" << int(t.first) << " " << debug::PrintComplexPosition(t.second));
        }
        DEBUG(debug::Debug(g_, rc_paths[i]));
        for (const auto &t : rc_threads[i]) {
            DEBUG("" << int(t.first) << " " << debug::PrintComplexPosition(t.second));
        }
    }
    {
        size_t i = chosen_path;
        DEBUG(debug::Debug(g_, paths[i]));
        for (const auto &t : threads_to_delete[i]) {
            DEBUG("" << int(t.first) << " " << debug::PrintComplexPosition(t.second));
        }
        DEBUG(debug::Debug(g_, rc_paths[i]));
        for (const auto &t : rc_threads[i]) {
            DEBUG("" << int(t.first) << " " << debug::PrintComplexPosition(t.second));
        }
    }

    LockDelete();
    coordinates_handler_.LockChanges();
    for (size_t i = 0; i < paths.size(); ++i) {
      /*
      if (threads_to_delete[i].size() != rc_threads[i].size()) {
        INFO("failed: paths " << Debug(paths[i]) << " and " << Debug(rc_paths[i]));
        for (const auto &e : threads_to_delete[i])
          INFO("to delete 1: " << int(e.first) << " - " << e.second);
        for (const auto &e : rc_threads[i])
          INFO("to delete 2: " << int(e.first) << " - " << e.second);
        INFO("---------paths1:");
        for (const auto &e : paths[i])
          coordinates_handler_.DebugOutput(e);
        INFO("---------paths2:");
        for (const auto &e : rc_paths[i])
          coordinates_handler_.DebugOutput(e);
        VERIFY(false);
      }
      */

      if (i == chosen_path) continue;
      if (bad_paths.count(i) > 0) continue;

      if (threads_to_delete[i].size() == 0) {
        TRACE("nothin to delete!");
        continue;
      }
      bool success = true;

      success &= ProjectPath(paths[i], paths[chosen_path], threads_to_delete[i]);
      success &= ProjectPath(rc_paths[i], rc_paths[chosen_path], rc_threads[i]);

      if (!success) {
          ClearDeleteList();
          coordinates_handler_.UnrollChanges();
          break;
      }
    }
    coordinates_handler_.ReleaseChanges();
    ReleaseDelete();

    TRACE("CollapsePaths End");
    return true;
  }

  bool ProjectPath(const Path &from, const Path &to,
      const PosArray &threads_to_delete) {
    const bool success = coordinates_handler_.ProjectPath(
            from, to, threads_to_delete);

    if (!success)
        return false;

    for (const auto e : from) {
      if (coordinates_handler_.GetMultiplicity(e) == 0) {
        DeleteEdge(e);
      }
    }

    return true;
  }

  bool ProjectPath(const Path &from, const Path &to) {
    const std::vector<std::pair<uint, size_t> > threads_to_delete =
        coordinates_handler_.GetContiguousThreads(from);
    return ProjectPath(from, to, threads_to_delete);
  }


 private:
  class DeletingMergeHandler {
   public:
    DeletingMergeHandler(const std::vector<EdgeId> &to_delete)
        : delete_set_(to_delete.begin(), to_delete.end()) {
    }

    bool operator()(const std::vector<EdgeId> &old_edges, EdgeId new_edge) {
      return false;
      bool ret = true;
      for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
        ret &= bool(delete_set_.count(*it));
        delete_set_.erase(*it);
      }
      if (ret) {
        delete_set_.insert(new_edge);
      }
      return ret;
    }

   private:
    std::unordered_set<EdgeId> delete_set_;
  };

  bool CheckForRCAbsence(
      const std::vector<std::vector<std::pair<uint, size_t> > > &threads) {
    return true;
    // We use that #of thread = 2 * genome_num + RC
    std::vector<uint> genomes;
    for (const auto &entry : threads) {
      for (const auto &pos : entry) {
        genomes.push_back(pos.first);
      }
    }
    std::sort(genomes.begin(), genomes.end());
    for (size_t i = 1; i < genomes.size(); ++i) {
      if (genomes[i] == genomes[i - 1] + 1 &&
          (genomes[i] & 1) == 1) {
        return false;
      }
    }

    return true;
  }

  bool CheckForCorrectPaths(const std::vector<Path> &paths) {
    for (size_t i = 1; i < paths.size(); ++i) {
      if (!coordinates_handler_.CheckCorrectPathProjection(paths[i - 1], 
                                                           paths[i]))
        return false;
    }
    return true;
  }

  std::unordered_set<size_t> GetIntersectingPaths(const std::vector<Path> &paths,
      const size_t chosen_path) {
    std::unordered_set<EdgeId> intouchable_edges;
    for (const auto e : paths[chosen_path]) {
      intouchable_edges.insert(e);
    }

    std::unordered_set<size_t> result;
    for (size_t i = 0; i < paths.size(); ++i) {
      if (i == chosen_path) continue;

      for (const auto e : paths[i]) {
        if (intouchable_edges.count(e) > 0) {
          result.insert(i);
          break;
        }
      }
    }

    return result;
  }

  bool CheckDeletionOfIntouchables(const std::vector<Path> &paths,
      const std::vector<std::vector<std::pair<uint, size_t> > > &del_threads,
      const size_t chosen_path) const {
    std::unordered_set<EdgeId> intouchable_edges;
    for (const auto e : paths[chosen_path]) {
      intouchable_edges.insert(e);
      intouchable_edges.insert(g_.conjugate(e));
    }
    
    for (size_t i = 0; i < paths.size(); ++i) {
      if (i == chosen_path) continue;
      const Path &path = paths[i];

      for (const auto e : path) {
        if (coordinates_handler_.GetMultiplicity(e) == del_threads[i].size()) {
          if (intouchable_edges.count(e) > 0) {
            return false;
          }
        }
      }
    }

    return true;
  }

  bool CheckPresenceOfSelfContiguousEdges(const std::vector<Path> &paths) const {
      for (const auto &path : paths) {
          for (const auto e : path) {
              if (g_.conjugate(e) == e)
                  return true;
          }
      }

      return false;
  }

  size_t CalcBridges(const Path &path, const size_t thin_multiplicity) const {
    size_t bridges = 0;
    for (const auto e : path)
      bridges += coordinates_handler_.GetMultiplicity(e) == thin_multiplicity;

    return bridges;
  }

  size_t CalcMultiplicitySum(const Path &path) const {
    size_t multiplicity_sum = 0;
    for (const auto e : path)
      multiplicity_sum += coordinates_handler_.GetMultiplicity(e);

    return multiplicity_sum;
  }

  void DeleteEdge(const EdgeId edge) {
    edges_to_delete_.push_back(edge);
    if (!is_deleting_locked_)
      ReleaseDelete();
  }

  void ForceDeleteEdges(const std::vector<EdgeId> &edges) {
    VERIFY(!is_deleting_locked_);

    //TRACE("DeleteEdges Begin " << Debug(edges) << " of size " << edges.size());
    DeletingMergeHandler merge_handler(edges);
    ConditionedSmartSetIterator<Graph, EdgeId, DeletingMergeHandler> smart_it(
        g_, edges.begin(), edges.end(), merge_handler);

    for (; !smart_it.IsEnd(); ++smart_it) {
      edge_remover_.DeleteEdge(*smart_it);
    }
    TRACE("DeleteEdges End");
  }

  void LockDelete() {
    is_deleting_locked_ = true;
  }
  void ReleaseDelete() {
    is_deleting_locked_ = false;
    ForceDeleteEdges(edges_to_delete_);
    ClearDeleteList();
  }
  void ClearDeleteList() {
    edges_to_delete_.clear();
  }
  template<class T>
  std::string Debug(const std::vector<T> &p) {
    std::stringstream ss;
    for (const auto &x : p) {
      ss << g_.str(x) << ";";
    }
    return ss.str();
  }

  Graph &g_;
  CoordinatesHandler<Graph> &coordinates_handler_;
  EdgeRemover<Graph> edge_remover_;

  bool is_deleting_locked_;
  std::vector<EdgeId> edges_to_delete_;
  
  DECL_LOGGER("PathProjector")
    ;
};

}

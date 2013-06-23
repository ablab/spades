#include <iostream>
#include <vector>

#include "coordinates_handler.hpp"
//#include "coloring.hpp"
#include "omni/graph_processing_algorithm.hpp"

namespace cap {

template <class Graph>
class PathProjector {
 public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
  typedef std::vector<EdgeId> Path;
  typedef unsigned char uchar;

  PathProjector(Graph &g, CoordinatesHandler<Graph> &coordinates_handler)
      : g_(g),
        coordinates_handler_(coordinates_handler),
        edge_remover_(g),
        is_deleting_locked_(false) {
  }

  bool CollapsePaths(const std::vector<Path> &paths) {
    TRACE("CollapsePaths Begin");

    std::vector<std::vector<std::pair<uchar, size_t> > > threads_to_delete;
    std::vector<size_t> num_bridges;
    std::vector<size_t> sum_multiplicities;
    for (const auto &path : paths) {
      threads_to_delete.push_back(coordinates_handler_.GetContiguousThreads(path));
      num_bridges.push_back(CalcBridges(path, threads_to_delete.back().size()));
      sum_multiplicities.push_back(CalcMultiplicitySum(path));
    }
    // RC paths
    std::vector<Path> rc_paths;
    std::vector<std::vector<std::pair<uchar, size_t> > > rc_threads;
    for (const auto &path : paths) {
      Path rc_path;
      for (auto it = path.rbegin(); it != path.rend(); ++it) {
        rc_path.push_back(g_.conjugate(*it));
      }
      rc_paths.push_back(rc_path);
      rc_threads.push_back(coordinates_handler_.GetContiguousThreads(rc_path));
    }

    // We restrict merging RC and original genome threads
    if (!CheckForRCAbsence(threads_to_delete))
      return false;
    if (!CheckForCorrectPaths(paths))
      return false;


    size_t chosen_path = 0;
    for (size_t i = 1; i < paths.size(); ++i) {
      if (num_bridges[i] < num_bridges[chosen_path] ||
            (num_bridges[i] == num_bridges[chosen_path] &&
            sum_multiplicities[i] > sum_multiplicities[chosen_path])) {
        chosen_path = i;
      }
    }

    LockDelete();
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
      if (threads_to_delete[i].size() == 0) {
        TRACE("nothin to delete!");
        continue;
      }
      ProjectPath(paths[i], paths[chosen_path], threads_to_delete[i]);
      ProjectPath(rc_paths[i], rc_paths[chosen_path], rc_threads[i]);
    }
    ReleaseDelete();

    TRACE("CollapsePaths End");
    return true;
  }

  void ProjectPath(const Path &from, const Path &to,
      const std::vector<std::pair<uchar, size_t> > &threads_to_delete) {
    coordinates_handler_.ProjectPath(from, to, threads_to_delete);

    for (const auto e : from) {
      if (coordinates_handler_.GetMultiplicity(e) == 0) {
        DeleteEdge(e);
      }
    }
  }

  inline void ProjectPath(const Path &from, const Path &to) {
    const std::vector<std::pair<uchar, size_t> > threads_to_delete =
        coordinates_handler_.GetContiguousThreads(from);
    ProjectPath(from, to, threads_to_delete);
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
      const std::vector<std::vector<std::pair<uchar, size_t> > > &threads) {
    return true;
    // We use that #of thread = 2 * genome_num + RC
    std::vector<uchar> genomes;
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


#pragma once
#include "common/pipeline/graph_pack.hpp"
#include <common/modules/path_extend/extension_chooser.hpp>
#include <common/modules/alignment/long_read_mapper.hpp>
#include "common/barcode_index/cluster_storage_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"

namespace path_extend {
namespace validation {

class ContigTransitionStorage {
 public:
    typedef path_extend::transitions::Transition Transition;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
 private:
    std::unordered_set<Transition> transitions_;
    std::unordered_set<EdgeId> covered_edges_;

 public:
    typedef std::unordered_set<Transition>::const_iterator const_iterator;
    void InsertTransition(const EdgeId& first, const EdgeId& second) {
        VERIFY(IsEdgeCovered(first) and IsEdgeCovered(second));
        transitions_.insert({{first, second}});
    }

    void InsertEdge(const EdgeId& edge) {
        covered_edges_.insert(edge);
    }

    const_iterator begin() const {
        return transitions_.begin();
    }

    const_iterator end() const {
        return transitions_.end();
    }

    size_t size() const {
        return transitions_.size();
    }

    bool CheckTransition(const ScaffoldVertex &first, const ScaffoldVertex &second) const {
        auto is_covered = [this](const EdgeId& edge) {
          return IsEdgeCovered(edge);
        };
        boost::optional<EdgeId> first_covered = first.getLastEdgeWithPredicate(is_covered);
        boost::optional<EdgeId> second_covered = second.getFirstEdgeWithPredicate(is_covered);
        if (first_covered.is_initialized() and second_covered.is_initialized()) {
            return CheckTransition(first_covered.get(), second_covered.get());
        }
        return false;
    }

    bool CheckTransition(const EdgeId& first, const EdgeId& second) const {
        Transition t(first, second);
        return CheckTransition(t);
    }

    bool CheckTransition(const Transition& transition) const {
        return transitions_.find(transition) != transitions_.end();
    }

    bool IsEdgeCovered(const EdgeId& edge) const {
        return covered_edges_.find(edge) != covered_edges_.end();
    }

    bool IsEdgeCovered(const ScaffoldVertex& vertex) const {
        auto is_covered = [this](const EdgeId& edge) {
          return IsEdgeCovered(edge);
        };
        boost::optional<EdgeId> first_covered = vertex.getFirstEdgeWithPredicate(is_covered);
        return first_covered.is_initialized();
    }

    bool CheckPath(const vector<EdgeId>& path) const {
        for (auto it1 = path.begin(), it2 = std::next(it1); it2 != path.end(); ++it1, ++it2) {
            EdgeId first = *it1;
            EdgeId second = *it2;
            if (not CheckTransition(first, second)) {
                return false;
            }
        }
        return true;
    }

    unordered_set<EdgeId> GetCoveredEdges() const {
        return covered_edges_;
    }
};


struct EdgeWithMapping {
  EdgeId edge_;
  Range mapping_;

  EdgeWithMapping(const EdgeId &edge_, const Range &mapping_) : edge_(edge_), mapping_(mapping_) {}
};

struct NamedPath {
  omnigraph::MappingPath<EdgeId> mapping_path;
  const string name;

  NamedPath(const MappingPath<EdgeId>& mapping_path, const string& name) : mapping_path(mapping_path), name(name) {}
};

struct NamedSimplePath {
  vector<EdgeWithMapping> path_;
  const string name_;

  NamedSimplePath(const vector<EdgeWithMapping>& path_, const string& name_) : path_(path_), name_(name_) {}
};

class ContigPathBuilder {
    const debruijn_graph::conj_graph_pack& gp_;

 public:
    ContigPathBuilder(const debruijn_graph::conj_graph_pack& graph_pack) : gp_(graph_pack) {}

    vector<NamedSimplePath> GetContigPaths(const string& path_to_contigs) const;

    vector<vector<EdgeWithMapping>> StripNames(const vector<NamedSimplePath>& named_paths) const;

 protected:
    vector<NamedPath> GetRawPaths(const string &contig_path) const;

    string RemoveSpacesFromName(const string& name) const;

    vector<NamedSimplePath> FixMappingPaths(const vector<NamedPath>& contig_paths) const;
};

class ContigPathFilter {
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;

 public:
    ContigPathFilter(const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_)
        : unique_storage_(unique_storage_) {}

    vector<vector<EdgeWithMapping>> FilterPathsUsingUniqueStorage(const vector<vector<EdgeWithMapping>>& paths) const;

    vector<vector<EdgeWithMapping>> FilterPathsUsingLength(const vector<vector<EdgeWithMapping>>& paths,
                                                           const size_t min_length, const Graph& g) const;

    vector<vector<EdgeWithMapping>> MergeSameEdges(const vector<vector<EdgeWithMapping>>& paths) const;
};

class FilteredReferencePathHelper {
 private:
    const debruijn_graph::conj_graph_pack& gp_;
 public:
    explicit FilteredReferencePathHelper(const conj_graph_pack& gp_);

    vector<vector<EdgeWithMapping>> GetFilteredReferencePathsFromLength(const string& path_to_reference, size_t length_threshold);
};

class TransitionStorageBuilder {

 public:
    ContigTransitionStorage GetTransitionStorage(const vector<vector<EdgeWithMapping>>& contig_paths) const {
        return BuildStorage(contig_paths);
    }

    virtual ~TransitionStorageBuilder() = default;

 protected:

    virtual ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_paths) const = 0;

    DECL_LOGGER("TransitionStorageBuilder");
};


class StrictTransitionStorageBuilder : public TransitionStorageBuilder {

 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override;
};

class ReverseTransitionStorageBuilder : public TransitionStorageBuilder {

 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override;
};

class ConjugateTransitionStorageBuilder: public TransitionStorageBuilder {
 protected:
    const Graph& g_;
 public:
    ConjugateTransitionStorageBuilder(const Graph& g_) : g_(g_) {}
 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override;
};

class GeneralTransitionStorageBuilder: public TransitionStorageBuilder {
 protected:
    const Graph& g_;
    size_t distance_;
    const bool with_reverse_;
    const bool with_conjugate_;

 public:
    GeneralTransitionStorageBuilder(const Graph& g_,
                                    size_t distance_,
                                    const bool with_reverse_,
                                    const bool with_conjugate_)
        : g_(g_), distance_(distance_), with_reverse_(with_reverse_), with_conjugate_(with_conjugate_) {}
 public:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& paths) const override;
 private:

    void ProcessPair(EdgeId first, EdgeId second, bool with_conjugate, bool with_reverse, ContigTransitionStorage& storage) const;
};

class ApproximateTransitionStorageBuilder : public TransitionStorageBuilder {

 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override;
};

class ClusterTransitionExtractor {
 public:
    typedef path_extend::transitions::Transition Transition;
 private:
    const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer_;

 public:
    explicit ClusterTransitionExtractor(const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer_) : ordering_analyzer_(
        ordering_analyzer_) {}

    vector<Transition> ExtractAllTransitionsFromNonPathCluster(const cluster_storage::Cluster& cluster);

    vector<Transition> ExtractGoodTransitionsFromNonPathCluster(const cluster_storage::Cluster& cluster);

    vector<Transition> ExtractTransitionsFromOrdering(const vector<EdgeId>& ordering);

    vector<Transition> ExtractTransitionsFromPathCluster(const cluster_storage::Cluster& cluster) {
        auto ordering = ordering_analyzer_.GetOrderingFromCluster(cluster);
        return ExtractTransitionsFromOrdering(ordering);
    }
};

}
}

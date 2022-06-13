//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/extension_chooser.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"
#include "alignment/edge_index.hpp"
#include "alignment/kmer_mapper.hpp"
#include "alignment/long_read_mapper.hpp"
#include "alignment/sequence_mapper.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

//fixme duplication with long_read_mapper.cpp. Discuss
class GappedPathExtractor {
  public:
    typedef debruijn_graph::Graph Graph;

    GappedPathExtractor(const Graph& g) noexcept;
    virtual ~GappedPathExtractor() = default;

    std::vector<PathWithMappingInfo> operator() (const MappingPath<EdgeId>& mapping) const;

  protected:

    virtual MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                  const MappingPath<EdgeId> &mapping_path) const;

    std::vector<PathWithMappingInfo> FindReadPathWithGaps(MappingPath<EdgeId> &path) const;

    size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path,
                               size_t& mapping_index, MappingRange& range_of_mapped_edge) const;

    const Graph& g_;
    const debruijn_graph::MappingPathFixer<Graph> path_fixer_;
    constexpr static double MIN_MAPPED_RATIO = 0.3;
    constexpr static size_t MIN_MAPPED_LENGTH = 100;
};

class ContigTransitionStorage {
    //todo Replace with scaffold graph?
  public:
    typedef transitions::Transition Transition;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_set<Transition>::const_iterator const_iterator;

    explicit ContigTransitionStorage(const ScaffoldingUniqueEdgeStorage &unique_storage) :
        transitions_(),
        covered_edges_(),
        unique_storage_(unique_storage) {}

    void InsertTransition(const EdgeId &first, const EdgeId &second) {
        VERIFY(IsEdgeCovered(first) and IsEdgeCovered(second));
        transitions_.insert({{first, second}});
    }
    void InsertEdge(const EdgeId &edge) {
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
        auto is_covered = [this](const EdgeId &edge) {
          return IsEdgeCovered(edge);
        };
        boost::optional<EdgeId> first_covered = first.GetLastEdgeWithPredicate(is_covered);
        boost::optional<EdgeId> second_covered = second.GetFirstEdgeWithPredicate(is_covered);
        if (first_covered.is_initialized() and second_covered.is_initialized()) {
            return CheckTransition(first_covered.get(), second_covered.get());
        }
        return false;
    }
    bool CheckTransition(const EdgeId &first, const EdgeId &second) const {
        Transition t(first, second);
        return CheckTransition(t);
    }
    bool CheckTransition(const Transition &transition) const {
        return transitions_.find(transition) != transitions_.end();
    }
    bool IsEdgeCovered(const EdgeId &edge) const {
        return covered_edges_.find(edge) != covered_edges_.end();
    }
    bool IsEdgeCovered(const ScaffoldVertex &vertex) const {
        auto is_covered = [this](const EdgeId &edge) {
          return IsEdgeCovered(edge);
        };
        boost::optional<EdgeId> first_covered = vertex.GetFirstEdgeWithPredicate(is_covered);
        return first_covered.is_initialized();
    }
    bool IsUnique(const EdgeId &edge) const {
        return unique_storage_.IsUnique(edge);
    }
    bool CheckPath(const std::vector<EdgeId> &path) const {
        for (auto it1 = path.begin(), it2 = std::next(it1); it2 != path.end(); ++it1, ++it2) {
            EdgeId first = *it1;
            EdgeId second = *it2;
            if (not CheckTransition(first, second)) {
                return false;
            }
        }
        return true;
    }
    std::unordered_set<EdgeId> GetCoveredEdges() const {
        return covered_edges_;
    }

  private:
    std::unordered_set<Transition> transitions_;
    std::unordered_set<EdgeId> covered_edges_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
};

struct EdgeWithMapping {
  EdgeWithMapping(const EdgeId &edge, const Range &mapping) : edge_(edge), mapping_(mapping) {}

  EdgeId edge_;
  Range mapping_;
};

struct NamedPath {
  NamedPath(const MappingPath<EdgeId> &mapping_path, const string &name) : mapping_path(mapping_path), name(name) {}

  omnigraph::MappingPath<EdgeId> mapping_path;
  const string name;

};

struct NamedSimplePath {
  NamedSimplePath(const std::vector<EdgeWithMapping> &path_, const string &name_) : path_(path_), name_(name_) {}

  std::vector<EdgeWithMapping> path_;
  const string name_;
};

class ContigPathBuilder {
  public:
    ContigPathBuilder(const Graph &g,
                      const debruijn_graph::EdgeIndex<Graph> &index,
                      const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                      size_t path_length_threshold) :
        g_(g),
        index_(index),
        kmer_mapper_(kmer_mapper),
        path_length_threshold_(path_length_threshold) {}

    std::vector<NamedSimplePath> GetContigPaths(const std::filesystem::path &path_to_contigs) const;

    std::vector<NamedPath> GetRawPaths(const std::filesystem::path &contig_path) const;

    std::vector<std::vector<EdgeWithMapping>> StripNames(const std::vector<NamedSimplePath> &named_paths) const;

  protected:
    string RemoveSpacesFromName(const string &name) const;

    std::vector<NamedSimplePath> FixMappingPaths(const std::vector<NamedPath> &contig_paths) const;

  private:
    const Graph &g_;
    const debruijn_graph::EdgeIndex<Graph> &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    size_t path_length_threshold_;
};

class ContigPathFilter {
  public:
    typedef std::vector<std::vector<EdgeWithMapping>> ReferencePaths;
    ContigPathFilter(const func::TypedPredicate<const EdgeId&> pred): pred_(pred) {}
    ReferencePaths FilterPaths(const ReferencePaths &paths) const;

  private:
    ReferencePaths MergeSameEdges(const ReferencePaths &paths) const;
    ReferencePaths RemoveRepeats(const ReferencePaths &paths) const;
    std::unordered_set<EdgeId> GetRepeats(const ReferencePaths &paths) const;

    const func::TypedPredicate<EdgeId> pred_;
};

struct UniqueReferencePaths {
  typedef std::vector<std::vector<EdgeWithMapping>> Paths;
  UniqueReferencePaths(Paths &&paths, const ScaffoldingUniqueEdgeStorage &unique_storage):
  paths_(paths), unique_storage_(unique_storage) {}

  size_t size() const {
      return paths_.size();
  }

  std::vector<std::vector<EdgeWithMapping>> paths_;
  const ScaffoldingUniqueEdgeStorage &unique_storage_;
};

class FilteredReferencePathHelper {
  public:
    typedef std::vector<std::vector<EdgeWithMapping>> ReferencePaths;
    FilteredReferencePathHelper(const Graph &g,
                                const debruijn_graph::EdgeIndex<Graph> &index,
                                const debruijn_graph::KmerMapper<Graph> &kmer_mapper) :
        g_(g),
        index_(index),
        kmer_mapper_(kmer_mapper) {};

    UniqueReferencePaths GetFilteredReferencePathsFromUnique(const std::filesystem::path &path_to_reference,
                                                             const ScaffoldingUniqueEdgeStorage &unique_storage) const;

  private:
    const Graph &g_;
    const debruijn_graph::EdgeIndex<Graph> &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
};

class TransitionStorageBuilder {

  public:
    ContigTransitionStorage GetTransitionStorage(const UniqueReferencePaths &contig_paths) const {
        return BuildStorage(contig_paths);
    }

    virtual ~TransitionStorageBuilder() = default;

  protected:

    virtual ContigTransitionStorage BuildStorage(const UniqueReferencePaths &contig_paths) const = 0;

    DECL_LOGGER("TransitionStorageBuilder");
};

class StrictTransitionStorageBuilder : public TransitionStorageBuilder {
  public:

  protected:
    ContigTransitionStorage BuildStorage(const UniqueReferencePaths &contig_paths) const override;
};

class ReverseTransitionStorageBuilder : public TransitionStorageBuilder {
  protected:
    ContigTransitionStorage BuildStorage(const UniqueReferencePaths &contig_paths) const override;
};

class ConjugateTransitionStorageBuilder : public TransitionStorageBuilder {
  public:
    ConjugateTransitionStorageBuilder(const Graph &g) : g_(g) {}
  protected:
    ContigTransitionStorage BuildStorage(const UniqueReferencePaths &contig_paths) const override;

    const Graph &g_;
};

class GeneralTransitionStorageBuilder : public TransitionStorageBuilder {
  public:
    GeneralTransitionStorageBuilder(const Graph &g_,
                                    size_t distance_,
                                    const bool with_reverse_,
                                    const bool with_conjugate_)
        : g_(g_), distance_(distance_), with_reverse_(with_reverse_), with_conjugate_(with_conjugate_) {}

    ContigTransitionStorage BuildStorage(const UniqueReferencePaths &contig_paths) const override;

  protected:
    const Graph &g_;
    size_t distance_;
    const bool with_reverse_;
    const bool with_conjugate_;

  private:
    void ProcessPair(EdgeId first,
                     EdgeId second,
                     bool with_conjugate,
                     bool with_reverse,
                     ContigTransitionStorage &storage) const;
};

class ClusterTransitionExtractor {
  public:
    typedef path_extend::read_cloud::transitions::Transition Transition;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    explicit ClusterTransitionExtractor(const cluster_storage::ClusterGraphAnalyzer &ordering_analyzer_)
        : ordering_analyzer_(
        ordering_analyzer_) {}

    std::vector<Transition> ExtractAllTransitionsFromNonPathCluster(const cluster_storage::Cluster &cluster);
    std::vector<Transition> ExtractGoodTransitionsFromNonPathCluster(const cluster_storage::Cluster &cluster);

  private:
    const cluster_storage::ClusterGraphAnalyzer &ordering_analyzer_;
};

}
}
}
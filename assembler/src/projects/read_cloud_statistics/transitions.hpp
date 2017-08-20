#pragma once
#include <common/pipeline/graph_pack.hpp>
#include <common/modules/path_extend/extension_chooser.hpp>
#include <common/modules/alignment/long_read_mapper.hpp>
#include "common/barcode_index/cluster_storage_extractor.hpp"
//#include "barcode_index/barcode_index.hpp"

namespace transitions {
struct Transition {
    EdgeId first_;
    EdgeId second_;
 public:
    Transition(const EdgeId &first, const EdgeId &second) : first_(first), second_(second) {}
    bool operator== (const Transition& other) const {
        return first_ == other.first_ and second_ == other.second_;
    };

    bool operator <(const Transition& other) const {
        return first_.int_id() < other.first_.int_id() or (first_.int_id() == other.first_.int_id() and
            second_.int_id() < other.second_.int_id());
    }

    Transition& operator =(const Transition& other) = default;
};
}

namespace std {
template<>
struct hash<transitions::Transition> {
  size_t operator()(const transitions::Transition &transition) const {
      using std::hash;
      return hash<size_t>()(transition.first_.int_id() + transition.second_.int_id());
  }
};
}

namespace transitions {
class ContigTransitionStorage {
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

    size_t Size() const {
        return transitions_.size();
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
};


struct EdgeWithMapping {
  EdgeId edge_;
  Range mapping_;

  EdgeWithMapping(const EdgeId &edge_, const Range &mapping_) : edge_(edge_), mapping_(mapping_) {}
};

class ContigPathBuilder {
    const debruijn_graph::conj_graph_pack& gp_;

 public:
    ContigPathBuilder(const debruijn_graph::conj_graph_pack& graph_pack) : gp_(graph_pack) {}

    vector<vector<EdgeWithMapping>> GetContigPaths(const string& path_to_contigs) const {
        auto raw_contig_paths = GetRawPaths(path_to_contigs);
        INFO(raw_contig_paths.size() << " raw paths");
        auto fixed_paths = FixMappingPaths(raw_contig_paths);
        INFO(fixed_paths.size() << " fixed paths");
        return fixed_paths;
    }

 protected:
    vector<omnigraph::MappingPath<EdgeId>> GetRawPaths(const string &contig_path) const {
        vector<omnigraph::MappingPath<EdgeId>> reference_paths;
        const debruijn_graph::Index &index = gp_.index;
        const debruijn_graph::KmerMapper<Graph>& kmer_mapper = gp_.kmer_mapper;
        auto contig_stream_ptr = make_shared<io::FileReadStream>(contig_path);
        auto rc_contig_stream_ptr = io::RCWrap<io::SingleRead>(contig_stream_ptr);

        vector<io::SingleRead> references;
        const size_t contig_length_threshold = 4000;
        while(!rc_contig_stream_ptr->eof()) {
            io::SingleRead contig;
            (*rc_contig_stream_ptr) >> contig;
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::Index>>
                (gp_.g, index, kmer_mapper);
            MappingPath<EdgeId> path = mapper->MapRead(contig);
            if (path.size() > 0 and path.back().second.initial_range.end_pos > contig_length_threshold) {
                reference_paths.push_back(path);
            }
        }
        return reference_paths;
    }

    vector<vector<EdgeWithMapping>> FixMappingPaths(const vector<omnigraph::MappingPath<EdgeId>>& reference_paths) const {
        vector<vector<EdgeWithMapping>> result;
        debruijn_graph::GappedPathExtractor gapped_path_extractor(gp_.g);
        for (const auto& path: reference_paths) {
            vector<vector<EdgeId>> fixed_paths = gapped_path_extractor(path);
            vector<EdgeWithMapping> long_path;
            size_t prefix_len = 0;
            for (const auto& fixed_path: fixed_paths) {
                for (const auto& edge: fixed_path) {
                    Range mapping(prefix_len, prefix_len + gp_.g.length(edge));
                    EdgeWithMapping ewm(edge, mapping);
                    long_path.push_back(ewm);
                    prefix_len += gp_.g.length(edge);
                }
            }
            result.push_back(long_path);
        }
        return result;
    }
};

class ContigPathFilter {
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;

 public:
    ContigPathFilter(const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_)
        : unique_storage_(unique_storage_) {}

    vector<vector<EdgeWithMapping>> FilterPathsUsingUniqueStorage(const vector<vector<EdgeWithMapping>>& paths) const {
        vector<vector<EdgeWithMapping>> result;
        for (const auto& path: paths) {
            vector<EdgeWithMapping> filtered_path;
            for (const auto& ewm: path) {
                if (unique_storage_.IsUnique(ewm.edge_)) {
                    filtered_path.push_back(ewm);
                }
            }
            if (not filtered_path.empty()) {
                result.push_back(filtered_path);
            }
        }
        return MergeSameEdges(result);
    }

    vector<vector<EdgeWithMapping>> FilterPathsUsingLength(const vector<vector<EdgeWithMapping>>& paths,
                                                           const size_t min_length, const Graph& g) const {
        vector<vector<EdgeWithMapping>> result;
        for (const auto& path: paths) {
            vector<EdgeWithMapping> filtered_path;
            for (const auto& ewm: path) {
                if (g.length(ewm.edge_) >= min_length) {
                    filtered_path.push_back(ewm);
                }
            }
            if (not filtered_path.empty()) {
                result.push_back(filtered_path);
            }
        }
        return MergeSameEdges(result);
    }

    vector<vector<EdgeWithMapping>> MergeSameEdges(const vector<vector<EdgeWithMapping>>& paths) const {
        vector<vector<EdgeWithMapping>> result;
        for (const auto& path: paths) {
            vector<EdgeWithMapping> current_path;
            for (const auto& ewm: path) {
                if (current_path.empty() or current_path.back().edge_ != ewm.edge_) {
                    current_path.push_back(ewm);
                } else {
                    if (not current_path.empty()) {
                        Range prev_mapping = current_path.back().mapping_;
                        Range new_mapping(prev_mapping.start_pos, ewm.mapping_.end_pos);
                        EdgeWithMapping new_ewm(ewm.edge_, new_mapping);
                        current_path.back() = new_ewm;
                    }
                }
            }

            //merge last and first edge for cycle paths
            auto last_ewm = current_path.back();
            auto first_ewm = current_path[0];
            if (last_ewm.edge_ == first_ewm.edge_ and current_path.size() > 1) {
                Range new_mapping(std::min(last_ewm.mapping_.start_pos, first_ewm.mapping_.start_pos),
                                  std::max(last_ewm.mapping_.end_pos, first_ewm.mapping_.end_pos));
                current_path[0].mapping_ = new_mapping;
                current_path.pop_back();
            }
            if (not current_path.empty()) {
                result.push_back(current_path);
            }
        }
        return result;
    }
};

class TransitionStorageBuilder {
//    using namespace debruijn_graph;

    const ContigPathBuilder& contig_path_builder_;
 public:
    TransitionStorageBuilder(const ContigPathBuilder& contig_path_builder) : contig_path_builder_(contig_path_builder) {}

    ContigTransitionStorage GetTransitionStorage(const vector<vector<EdgeWithMapping>>& contig_paths) const {

        INFO(contig_paths.size() << " paths");

        return BuildStorage(contig_paths);
    }

    virtual ~TransitionStorageBuilder() = default;

 protected:

    virtual ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_paths) const = 0;

    DECL_LOGGER("TransitionStorageBuilder");
};

class StrictTransitionStorageBuilder : public TransitionStorageBuilder {
 public:
    StrictTransitionStorageBuilder(const ContigPathBuilder& contig_path_builder)
        : TransitionStorageBuilder(contig_path_builder) {}

 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override {
        ContigTransitionStorage storage;
        for (const auto& path: long_edges) {
            for (auto it1 = path.begin(), it2 = std::next(path.begin());
                 it1 != path.end() and it2 != path.end(); ++it1, ++it2) {
                EdgeId edge1 = (*it1).edge_;
                EdgeId edge2 = (*it2).edge_;
                storage.InsertEdge(edge1);
                storage.InsertEdge(edge2);
                storage.InsertTransition(edge1, edge2);
            }
        }
        return storage;
    }
};

class ApproximateTransitionStorageBuilder : public TransitionStorageBuilder {
 public:
    ApproximateTransitionStorageBuilder(const ContigPathBuilder& contig_path_builder)
        : TransitionStorageBuilder(contig_path_builder) {}

 protected:
    ContigTransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override {
        ContigTransitionStorage storage;
        const size_t threshold = 5000;
        const size_t mapping_error_threshold = 500;
        for (const auto& path: long_edges) {
            for (const auto& ewm: path) {
                storage.InsertEdge(ewm.edge_);
            }
        }
        for (const auto& path: long_edges) {
            for (auto it = path.begin(); it != path.end(); ++it) {
                EdgeId first = (*it).edge_;
                size_t pos = (*it).mapping_.end_pos;
                for (auto it_next = std::next(it); it_next != path.end(); ++it_next) {
                    EdgeId second = (*it_next).edge_;
                    size_t next_pos = (*it_next).mapping_.start_pos;
                    VERIFY(next_pos + mapping_error_threshold >= pos);
                    if (pos + threshold > next_pos) {
                        storage.InsertTransition(first, second);
                    } else {
                        break;
                    }
                }
            }
        }
        return storage;
    }
};

class ClusterTransitionExtractor {
    const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer_;

 public:
    ClusterTransitionExtractor(const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer_) : ordering_analyzer_(
        ordering_analyzer_) {}

    vector<transitions::Transition> ExtractTransitionsFromNonPathCluster(const cluster_storage::Cluster& cluster) {
        const auto& internal_graph = cluster.GetInternalGraph();
        vector<transitions::Transition> result;
        for (const auto& vertex: internal_graph) {
            for (auto it = internal_graph.outcoming_begin(vertex); it != internal_graph.outcoming_end(vertex); ++it) {
                result.emplace_back(vertex, *it);
            }
        }
        return result;
    }

    vector<transitions::Transition> ExtractTransitionsFromOrdering(const vector<EdgeId>& ordering) {
        VERIFY(ordering.size() != 0);
        vector<transitions::Transition> result;
        for (auto first = ordering.begin(), second = std::next(ordering.begin()); second != ordering.end(); ++first, ++second) {
            result.emplace_back(*first, *second);
        }
        return result;
    }

    vector<transitions::Transition> ExtractTransitionsFromPathCluster(const cluster_storage::Cluster& cluster) {
        auto ordering = ordering_analyzer_.GetOrderingFromCluster(cluster);
        return ExtractTransitionsFromOrdering(ordering);
    }
};
}




#pragma once
#include <common/pipeline/graph_pack.hpp>
#include <common/modules/path_extend/extension_chooser.hpp>
//#include "barcode_index/barcode_index.hpp"

namespace transitions {
struct Transition {
    const EdgeId first_;
    const EdgeId second_;
 public:
    Transition(const EdgeId &first, const EdgeId &second) : first_(first), second_(second) {}
    bool operator== (const Transition& other) const {
        return first_ == other.first_ and second_ == other.second_;
    };
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
class TransitionStorage {
    std::unordered_set<Transition> transitions_;

 public:
    typedef std::unordered_set<Transition>::const_iterator const_iterator;
    void Insert(const EdgeId& first, const EdgeId& second) {
        transitions_.insert({{first, second}});
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
        return (transitions_.find(t) != transitions_.end());
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

class CorrectTransitionsStorageBuilder {
//    using namespace debruijn_graph;

    const debruijn_graph::conj_graph_pack& gp_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
 public:
    CorrectTransitionsStorageBuilder(const debruijn_graph::conj_graph_pack& graph_pack,
                              const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage) :
        gp_(graph_pack), unique_storage_(unique_storage) {}

    TransitionStorage GetTransitionStorage(const string& genome_path) {
        INFO("Mapping references")
        auto reference_paths = GetReferencePaths(genome_path);
        INFO("Building transition storage");
        auto long_paths = GetLongPaths(reference_paths);
        return BuildStorage(long_paths);
    }

 protected:
    vector<vector<EdgeWithMapping>> GetLongPaths (const vector<omnigraph::MappingPath<EdgeId>>& reference_paths) const {
        vector<vector<EdgeWithMapping>> long_paths;
        for (const auto& path: reference_paths) {
            auto simple_path = path.simple_path();
            vector<EdgeWithMapping> path_of_long;
            for (size_t i = 0; i < path.size(); ++i) {
                EdgeWithMapping edge(path[i].first, path[i].second.initial_range);
                if (unique_storage_.IsUnique(edge.edge_)) {
                    path_of_long.push_back(edge);
                }
            }
            long_paths.push_back(path_of_long);
        }
        return long_paths;
    }

    virtual TransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_paths) const = 0;

    vector<omnigraph::MappingPath<EdgeId>> GetReferencePaths(const string& genome_path) const {
        vector<omnigraph::MappingPath<EdgeId>> reference_paths;
        const debruijn_graph::Index &index = gp_.index;
        const debruijn_graph::KmerMapper<Graph>& kmer_mapper = gp_.kmer_mapper;
        io::FileReadStream genome_stream(genome_path);
        vector<io::SingleRead> references;
        DEBUG(genome_path);
        while(!genome_stream.eof()) {
            io::SingleRead reference;
            DEBUG("reading reference");
            genome_stream >> reference;
            auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::Index>>
                (gp_.g, index, kmer_mapper);
            DEBUG("Mapping read");
            MappingPath<EdgeId> path = mapper->MapRead(reference);
            DEBUG("Size: " << path.size());
            if (path.size() > 0) {
                DEBUG(path.back().first.int_id());
                DEBUG(path.back().second.initial_range);
            }
            if (path.size() > 0 and path.back().second.initial_range.end_pos > 50000) {
                DEBUG("pushing path back");
                reference_paths.push_back(path);
            }
        }
        DEBUG("Returning");
        return reference_paths;
    }

    DECL_LOGGER("TransitionStorageBuilder");
};

class StrictTransitionStorageBuilder : public CorrectTransitionsStorageBuilder {
 public:
    StrictTransitionStorageBuilder(const debruijn_graph::conj_graph_pack &graph_pack,
                                   const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage)
        : CorrectTransitionsStorageBuilder(graph_pack, unique_storage) {

    }

 protected:
    TransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override {
        TransitionStorage storage;
        for (const auto& path: long_edges) {
            for (auto it1 = path.begin(), it2 = std::next(path.begin());
                 it1 != path.end() and it2 != path.end(); ++it1, ++it2) {
//                INFO("Iterators");
                EdgeId edge1 = (*it1).edge_;
                EdgeId edge2 = (*it2).edge_;
//                INFO("Inserting");
                storage.Insert(edge1, edge2);
            }
        }
        return storage;
    }
};

class ApproximateTransitionStorageBuilder : public CorrectTransitionsStorageBuilder {
 public:
    ApproximateTransitionStorageBuilder(const debruijn_graph::conj_graph_pack &graph_pack,
                                   const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage)
        : CorrectTransitionsStorageBuilder(graph_pack, unique_storage) {}

 protected:
    TransitionStorage BuildStorage(const vector<vector<EdgeWithMapping>>& long_edges) const override {
        TransitionStorage storage;
        const size_t threshold = 10000;
        const size_t mapping_error_threshold = 500;
        for (const auto& path: long_edges) {
            for (auto it = path.begin(); it != path.end(); ++it) {
                EdgeId first = (*it).edge_;
                size_t pos = (*it).mapping_.end_pos;
//                INFO("previous end: " << pos);
                for (auto it_next = std::next(it); it_next != path.end(); ++it_next) {
                    EdgeId second = (*it_next).edge_;
                    size_t next_pos = (*it_next).mapping_.start_pos;
//                    INFO("current start: " << next_pos);
                    VERIFY(next_pos + mapping_error_threshold >= pos);
                    if (pos + threshold > next_pos) {
                        storage.Insert(first, second);
                    } else {
                        break;
                    }
                }
            }
        }
        return storage;
    }
};
}




//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_vertex.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_extractor.hpp"

namespace path_extend {
namespace read_cloud {
namespace transitions {
struct Transition {
  public:
    typedef debruijn_graph::EdgeId EdgeId;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    Transition(const EdgeId &first, const EdgeId &second) : first_(first), second_(second) {}

    Transition(const ScaffoldVertex &first, const ScaffoldVertex &second) :
        first_(first.GetLastEdge()), second_(second.GetFirstEdge()) {}

    bool operator==(const Transition &other) const {
        return first_ == other.first_ and second_ == other.second_;
    };

    bool operator<(const Transition &other) const {
        return first_.int_id() < other.first_.int_id() or (first_.int_id() == other.first_.int_id() and
            second_.int_id() < other.second_.int_id());
    }

    Transition &operator=(const Transition &other) = default;

    EdgeId first_;
    EdgeId second_;
};
}
}
}

namespace std {
template<>
struct hash<path_extend::read_cloud::transitions::Transition> {
  size_t operator()(const path_extend::read_cloud::transitions::Transition &transition) const {
      using std::hash;
      return hash<size_t>()(transition.first_.int_id() + transition.second_.int_id());
  }
};
}

namespace path_extend {
namespace read_cloud {
namespace transitions {
typedef std::unordered_map<Transition, size_t> ClusterTransitionStorage;
}
}
}
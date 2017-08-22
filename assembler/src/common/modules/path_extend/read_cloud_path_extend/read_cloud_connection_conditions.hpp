#pragma once
#include "common/modules/path_extend/scaffolder2015/connection_condition2015.hpp"

namespace path_extend {
    //Same as AssemblyGraphConnectionCondition, but stops after reaching unique edges.
    class AssemblyGraphUniqueConnectionCondition : public AssemblyGraphConnectionCondition {
        using AssemblyGraphConnectionCondition::g_;
        using AssemblyGraphConnectionCondition::interesting_edge_set_;
        using AssemblyGraphConnectionCondition::max_connection_length_;
        //fixme duplication with interesting edges, needed to pass to dijkstra
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
     public:
        AssemblyGraphUniqueConnectionCondition(const Graph& g,
                                               size_t max_connection_length,
                                               const ScaffoldingUniqueEdgeStorage& unique_edges);
        map<EdgeId, double> ConnectedWith(EdgeId e) const override;
    };
}
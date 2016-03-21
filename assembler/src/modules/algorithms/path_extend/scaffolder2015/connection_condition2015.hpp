
#ifndef CONNECTION_CONDITION2015_HPP
#define CONNECTION_CONDITION2015_HPP
#include "algorithms/genome_consistance_checker.hpp"
#include "dev_support/logger/logger.hpp"
#include "algorithms/path_extend/paired_library.hpp"
#include "data_structures/assembly_graph/graph_support/scaff_supplementary.hpp"
#include <map>
#include <set>

namespace path_extend {

/* Connection condition are used by both scaffolder's extension chooser and scaffold graph */

class ConnectionCondition {
public:
// Outputs the edges e is connected with.
//TODO  performance issue: think about inside filtering. Return only unique connected edges?
    virtual set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const = 0;
// Outputs the weight of the pair e1 and e2
    virtual double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const = 0;
    virtual size_t GetLibIndex() const = 0;
    virtual ~ConnectionCondition() {
    }
};

// Main (mate pair library) connection condition.
class PairedLibConnectionCondition : public ConnectionCondition {
protected:
    const debruijn_graph::Graph &graph_;
    shared_ptr <PairedInfoLibrary> lib_;
    size_t lib_index_;
//Minimal number of mate pairs to call connection sound
    size_t min_read_count_;
public:
//Only paired info with gap between e1 and e2 between -left_dist_delta_ and right_dist_delta_ taken in account
    int left_dist_delta_;
    int right_dist_delta_;

    PairedLibConnectionCondition(const debruijn_graph::Graph &graph,
                                 shared_ptr <PairedInfoLibrary> lib,
                                 size_t lib_index,
                                 size_t min_read_count);
    size_t GetLibIndex() const override;
    set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const override;
    double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const override;
//Returns median gap size
    int GetMedianGap (debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const;
};

//Advanced mate-pair connection condition
class AdvancedPairedConnectionCondition: public PairedLibConnectionCondition {
protected:
    size_t always_add_;
    size_t never_add_;
    double relative_threshold_;

public:
    AdvancedPairedConnectionCondition(const debruijn_graph::Graph &graph,
                                      shared_ptr <PairedInfoLibrary> lib,
                                      size_t lib_index,
                                      size_t always_add,
                                      size_t never_add,
                                      double relative_threshold);

    set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const override;

};

/*  Condition used to find connected in graph edges.
*
*/
class AssemblyGraphConnectionCondition : public ConnectionCondition {
protected:
    const debruijn_graph::Graph &g_;
//Maximal gap to the connection.
    size_t max_connection_length_;
    set<EdgeId> interesting_edge_set_;
    mutable map <debruijn_graph::Graph::EdgeId, set<debruijn_graph::Graph::EdgeId>> stored_distances_;
public:
    AssemblyGraphConnectionCondition(const debruijn_graph::Graph &g, size_t max_connection_length,
                                     const ScaffoldingUniqueEdgeStorage& unique_edges);
    void AddInterestingEdge(debruijn_graph::EdgeId e);
    set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const override;
    double GetWeight(debruijn_graph::EdgeId, debruijn_graph::EdgeId) const override;
    size_t GetLibIndex() const override;
};
}

#endif //PROJECT_CONNECTION_CONDITION2015_HPP

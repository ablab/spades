
#ifndef CONNECTION_CONDITION2015_HPP
#define CONNECTION_CONDITION2015_HPP
#include "algorithms/genome_consistance_checker.hpp"
#include "dev_support/logger/logger.hpp"
#include "algorithms/path_extend/paired_library.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "assembly_graph/graph_alignment/long_read_storage.hpp"
#include "algorithms/path_extend/pe_utils.hpp"
#include <map>
#include <set>


namespace path_extend {

/* Connection condition are used by both scaffolder's extension chooser and scaffold graph */

class ConnectionCondition {
protected:
    DECL_LOGGER("ConnectionCondition")

public:
// Outputs the edges e is connected with.
//TODO  performance issue: think about inside filtering. Return only unique connected edges?
    virtual map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e) const = 0;
    virtual map <debruijn_graph::EdgeId, double > ConnectedWith(debruijn_graph::EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const;
    virtual int GetMedianGap (debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const = 0;
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
    map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e) const override;
    double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const;
//Returns median gap size
    int GetMedianGap (debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const override;
};

class LongReadsLibConnectionCondition : public ConnectionCondition {
protected:
    const debruijn_graph::Graph &graph_;
    size_t lib_index_;
//Minimal number of reads to call connection sound
    size_t min_read_count_;
    const GraphCoverageMap& cov_map_;
public:
//Only paired info with gap between e1 and e2 between -left_dist_delta_ and right_dist_delta_ taken in account

    LongReadsLibConnectionCondition(const debruijn_graph::Graph &graph,
                                 size_t lib_index,
                                 size_t min_read_count, const GraphCoverageMap& cov_map);
    size_t GetLibIndex() const override;
    map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e) const override;
    map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const override;
// Returns median gap size
    int GetMedianGap (debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const override;
};



//Should it be removed after ConnectedWith using unique storage was introduced?
class ScaffoldGraphPairedConnectionCondition: public PairedLibConnectionCondition {
protected:
    const set<debruijn_graph::EdgeId>& graph_edges_;

    size_t always_add_;
    size_t never_add_;
    double relative_threshold_;

public:
    ScaffoldGraphPairedConnectionCondition(const debruijn_graph::Graph &graph,
                                      const set<debruijn_graph::EdgeId>& graph_edges,
                                      shared_ptr <PairedInfoLibrary> lib,
                                      size_t lib_index,
                                      size_t always_add,
                                      size_t never_add,
                                      double relative_threshold);

    map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e) const override;

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
    mutable map <debruijn_graph::Graph::EdgeId, map<debruijn_graph::Graph::EdgeId, double>> stored_distances_;
public:
    AssemblyGraphConnectionCondition(const debruijn_graph::Graph &g, size_t max_connection_length,
                                     const ScaffoldingUniqueEdgeStorage& unique_edges);
    void AddInterestingEdge(debruijn_graph::EdgeId e);
    map <debruijn_graph::EdgeId, double> ConnectedWith(debruijn_graph::EdgeId e) const;
    size_t GetLibIndex() const override;
    int GetMedianGap (debruijn_graph::EdgeId, debruijn_graph::EdgeId ) const override;
};
}

#endif //PROJECT_CONNECTION_CONDITION2015_HPP

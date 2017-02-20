#pragma once
#include "modules/genome_consistance_checker.hpp"
#include "utils/logger/logger.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "modules/path_extend/pe_utils.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include <map>
#include <set>


namespace path_extend {
using debruijn_graph::EdgeId;
using debruijn_graph::Graph;

//De Bruijn graph edge condition interface
class LengthLowerBound : public omnigraph::EdgeCondition<Graph> {
    typedef Graph::EdgeId EdgeId;
    typedef EdgeCondition<Graph> base;

    const size_t max_length_;

public:

    LengthLowerBound(const Graph &g, size_t max_length)
            : base(g),
              max_length_(max_length) {
    }

    bool Check(EdgeId e) const {
        return this->g().length(e) >= max_length_;
    }
};

/* Connection condition are used by both scaffolder's extension chooser and scaffold graph */

class ConnectionCondition {
protected:
    DECL_LOGGER("ConnectionCondition")

public:
// Outputs the edges e is connected with.
//TODO  performance issue: think about inside filtering. Return only unique connected edges?
    virtual map<EdgeId, double> ConnectedWith(EdgeId e) const = 0;
    virtual map<EdgeId, double> ConnectedWith(EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const;
    virtual int GetMedianGap(EdgeId e1, EdgeId e2) const = 0;
    virtual size_t GetLibIndex() const = 0;
    virtual ~ConnectionCondition() {
    }
};

// Main (mate pair library) connection condition.
class PairedLibConnectionCondition : public ConnectionCondition {
protected:
    const Graph &graph_;
    shared_ptr <PairedInfoLibrary> lib_;
    size_t lib_index_;
//Minimal number of mate pairs to call connection sound
    size_t min_read_count_;
public:
//Only paired info with gap between e1 and e2 between -left_dist_delta_ and right_dist_delta_ taken in account
    int left_dist_delta_;
    int right_dist_delta_;

    PairedLibConnectionCondition(const Graph &graph,
                                 shared_ptr <PairedInfoLibrary> lib,
                                 size_t lib_index,
                                 size_t min_read_count);
    size_t GetLibIndex() const override;
    map <EdgeId, double> ConnectedWith(EdgeId e) const override;
    double GetWeight(EdgeId e1, EdgeId e2) const;
//Returns median gap size
    int GetMedianGap (EdgeId e1, EdgeId e2) const override;
};

class LongReadsLibConnectionCondition : public ConnectionCondition {
protected:
    const Graph &graph_;
    size_t lib_index_;
//Minimal number of reads to call connection sound
    size_t min_read_count_;
    const GraphCoverageMap& cov_map_;

    bool CheckPath(BidirectionalPath *path, EdgeId e1, EdgeId e2) const;

public:
//Only paired info with gap between e1 and e2 between -left_dist_delta_ and right_dist_delta_ taken in account

    LongReadsLibConnectionCondition(const Graph &graph,
                                 size_t lib_index,
                                 size_t min_read_count, const GraphCoverageMap& cov_map);
    size_t GetLibIndex() const override;
    map<EdgeId, double> ConnectedWith(EdgeId e) const override;
    map<EdgeId, double> ConnectedWith(EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const override;
// Returns median gap size
    int GetMedianGap (EdgeId e1, EdgeId e2) const override;

};



//Should it be removed after ConnectedWith using unique storage was introduced?
class ScaffoldGraphPairedConnectionCondition: public PairedLibConnectionCondition {
protected:
    const set<EdgeId>& graph_edges_;

    size_t always_add_;
    size_t never_add_;
    double relative_threshold_;

public:
    ScaffoldGraphPairedConnectionCondition(const Graph &graph,
                                      const set<EdgeId>& graph_edges,
                                      shared_ptr <PairedInfoLibrary> lib,
                                      size_t lib_index,
                                      size_t always_add,
                                      size_t never_add,
                                      double relative_threshold);

    map<EdgeId, double> ConnectedWith(EdgeId e) const override;

};

/*  Condition used to find connected in graph edges.
*
*/
class AssemblyGraphConnectionCondition : public ConnectionCondition {
protected:
    const Graph &g_;
//Maximal gap to the connection.
    size_t max_connection_length_;
    set<EdgeId> interesting_edge_set_;
    mutable map<EdgeId, map<EdgeId, double>> stored_distances_;
public:
    AssemblyGraphConnectionCondition(const Graph &g, size_t max_connection_length,
                                     const ScaffoldingUniqueEdgeStorage& unique_edges);
    void AddInterestingEdges(func::TypedPredicate<typename Graph::EdgeId> edge_condition);
    map<EdgeId, double> ConnectedWith(EdgeId e) const override;
    size_t GetLibIndex() const override;
    int GetMedianGap(EdgeId, EdgeId ) const override;
};
}

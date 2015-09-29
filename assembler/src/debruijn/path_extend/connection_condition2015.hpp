//
// Created by lab42 on 9/29/15.
//

#ifndef PROJECT_CONNECTION_CONDITION2015_HPP
#define PROJECT_CONNECTION_CONDITION2015_HPP
#include "genome_consistance_checker.hpp"
#include "logger/logger.hpp"
#include "paired_library.hpp"
#include <map>
#include <set>

namespace path_extend {

    class ConnectionCondition {
    public:
        virtual set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const = 0;
        virtual double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const = 0;
        virtual size_t GetLibIndex() const = 0;
        virtual ~ConnectionCondition() {
        }
    };

    class PairedLibConnectionCondition : public ConnectionCondition {
    private:
        const debruijn_graph::Graph &graph_;
        shared_ptr <PairedInfoLibrary> lib_;
        size_t lib_index_;
        size_t min_read_count_;
        int left_dist_delta_;
        int right_dist_delta_;

    public:
        PairedLibConnectionCondition(const debruijn_graph::Graph &graph,
                                     shared_ptr <PairedInfoLibrary> lib,
                                     size_t lib_index,
                                     size_t min_read_count);
        size_t GetLibIndex() const;
        set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const;
        double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const;
    };

    class AssemblyGraphConnectionCondition : public ConnectionCondition {
    private:
        const debruijn_graph::Graph &g_;
        size_t max_connection_length_;

    public:
        AssemblyGraphConnectionCondition(const debruijn_graph::Graph &g, size_t max_connection_length);

        set <debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const;
        double GetWeight(debruijn_graph::EdgeId, debruijn_graph::EdgeId) const;
        size_t GetLibIndex() const;
    };
}

#endif //PROJECT_CONNECTION_CONDITION2015_HPP

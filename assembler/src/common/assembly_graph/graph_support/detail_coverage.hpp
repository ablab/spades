//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/core/action_handlers.hpp"
#include "utils/verify.hpp"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <fstream>

namespace omnigraph {

template<class Graph>
class FlankingCoverage : public omnigraph::GraphActionHandler<Graph> {

public:
    typedef typename Graph::EdgeId EdgeId;

private:
    typedef omnigraph::GraphActionHandler<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef std::pair<EdgeId, unsigned> Pos;

    Graph& g_;
    const size_t averaging_range_;

public:
    void SetRawCoverage(EdgeId e, unsigned cov) {
        g_.data(e).set_flanking_coverage(cov);
    }

    unsigned RawCoverage(EdgeId e) const {
        return g_.data(e).flanking_coverage();
    }

private:
    size_t EdgeAveragingRange(EdgeId e) const {
        return std::min(this->g().length(e), averaging_range_);
    }

    double AverageFlankingCoverage(EdgeId e) const {
        return double(RawCoverage(e)) / double(EdgeAveragingRange(e));
    }

    unsigned InterpolateCoverage(EdgeId e, size_t l) const {
        VERIFY(l <= averaging_range_);
        VERIFY(l < g_.length(e));
        return unsigned(math::round(AverageFlankingCoverage(e) * double(l)));
    }

    void SetCoverageSimilarToAverageFlanking(EdgeId target, EdgeId source) {
        SetRawCoverage(target, unsigned(math::round(AverageFlankingCoverage(source) * double(EdgeAveragingRange(target)))));
    }

    void SetCoverageSimilarToAverageGlobal(EdgeId target, EdgeId source) {
        SetRawCoverage(target, unsigned(math::round(g_.coverage(source) * double(EdgeAveragingRange(target)))));
    }

public:
    //todo think about interactions with gap closer
    FlankingCoverage(Graph& g, size_t averaging_range)
            : base(g, "FlankingCoverage"), g_(g),
              averaging_range_(averaging_range) {
    }

    size_t averaging_range() const {
        return averaging_range_;
    }

    void IncRawCoverage(EdgeId e, unsigned count) {
        g_.data(e).inc_flanking_coverage(count);
    }

    double CoverageOfStart(EdgeId e) const {
        return AverageFlankingCoverage(e);
    }

    double CoverageOfEnd(EdgeId e) const {
        return CoverageOfStart(this->g().conjugate(e));
    }

    virtual void HandleAdd(EdgeId /*e*/) {
    }

    virtual void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) {
//        SetRawCoverage(new_edge, RawCoverage(old_edges.front()));
        size_t kpomers_left = averaging_range_;
        unsigned acc = 0;
        for (EdgeId e : old_edges) {
            if (kpomers_left >= g_.length(e)) {
                acc += RawCoverage(e);
                kpomers_left -= g_.length(e);
            } else {
                if (kpomers_left != 0)
                    acc += InterpolateCoverage(e, kpomers_left);
                break;
            }
        }
        SetRawCoverage(new_edge, acc);
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
        SetRawCoverage(new_edge, RawCoverage(edge1) + RawCoverage(edge2));
    }

    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                             EdgeId new_edge_2) {
        //todo maybe improve later
        SetCoverageSimilarToAverageFlanking(new_edge_1, old_edge);
        SetCoverageSimilarToAverageGlobal(new_edge_2, old_edge);
        if (old_edge == g_.conjugate(old_edge)) {
            SetCoverageSimilarToAverageGlobal(g_.conjugate(new_edge_1), old_edge);
        }
    }

    virtual void HandleDelete(EdgeId e) {
        SetRawCoverage(e, 0);
    }

    double LocalCoverage(EdgeId e, VertexId v) const {
        if (this->g().EdgeStart(e) == v) {
            return GetInCov(e);
        } else if (this->g().EdgeEnd(e) == v) {
            return GetOutCov(e);
        } else {
            VERIFY(false);
            return 0.0;
        }
    }

    //left for compatibility
    //todo rename
    double GetInCov(EdgeId e) const {
        return CoverageOfStart(e);
    }

    //left for compatibility
    //todo rename
    double GetOutCov(EdgeId e) const {
        return CoverageOfEnd(e);
    }

    //////////////////////////

    void Save(EdgeId e, std::ostream &out) const {
        out << RawCoverage(e);
    }

    void Load(EdgeId e, std::istream &in) {
        unsigned cov;
        in >> cov;
        SetRawCoverage(e, cov);
    }

    /*
     * Is thread safe if different threads process different edges.
     */
    bool IsThreadSafe() const {
        return true;
    }

private:
    DECL_LOGGER("FlankingCoverage");
};

}

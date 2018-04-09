//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * coverage.hpp
 *
 *  Created on: Jun 21, 2011
 *      Author: sergey
 */

#pragma once

#include "utils/logger/logger.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include "math/xmath.h"
#include "action_handlers.hpp"

namespace omnigraph {

//todo save/load raw k-mer coverage
template<class Graph>
class CoverageIndex : public GraphActionHandler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    Graph& g_;

 public:
    CoverageIndex(Graph &g)
            : GraphActionHandler<Graph>(g, "CoverageIndex"), g_(g) {
    }

    /**
     * In NON averaged units
     */
    void SetRawCoverage(EdgeId e, unsigned cov) {
        g_.data(e).set_raw_coverage(cov);
    }

    void IncRawCoverage(EdgeId e, unsigned count) {
        g_.data(e).inc_raw_coverage((int)count);
    }

    void SetAvgCoverage(EdgeId e, double cov) {
        g_.data(e).set_raw_coverage((int) math::round(cov * (double) this->g().length(e)));
    }

    /**
     * Returns average coverage of the edge
     */
    double coverage(EdgeId edge) const {
        return (double) RawCoverage(edge) / (double) this->g().length(edge);
    }

    unsigned RawCoverage(EdgeId edge) const {
        return g_.data(edge).raw_coverage();
    }

    void HandleDelete(EdgeId edge) override {
        SetRawCoverage(edge, 0);
    }

    void HandleMerge(const std::vector<EdgeId>& old_edges, EdgeId new_edge) override {
        unsigned coverage = 0;
        for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
            coverage += RawCoverage(*it);
        }
        SetRawCoverage(new_edge, coverage);
    }

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override {
        SetRawCoverage(new_edge, RawCoverage(edge1) + RawCoverage(edge2));
    }

    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override {
        double avg_cov = coverage(old_edge);
        if (old_edge == g_.conjugate(old_edge)) {
            int raw1 = std::max(1, (int) math::round(avg_cov * (double) this->g().length(new_edge1)));
            SetRawCoverage(new_edge1, raw1);
            SetRawCoverage(g_.conjugate(new_edge1), raw1);
            SetRawCoverage(new_edge2, std::max(1, (int) math::round(avg_cov * (double) this->g().length(new_edge2))));
        } else {
            SetRawCoverage(new_edge1, std::max(1, (int) math::round(avg_cov * (double) this->g().length(new_edge1))));
            SetRawCoverage(new_edge2, std::max(1, (int) math::round(avg_cov * (double) this->g().length(new_edge2))));
        }
    }

    void Save(EdgeId e, std::ostream& out) const {
        out << fmt::format("{:.6f}", coverage(e));
    }

    void Load(EdgeId e, std::istream& in) {
        double cov;
        in >> cov;
        SetAvgCoverage(e, cov);
    }

    /*
     * Is thread safe if different threads process different edges.
     */
    bool IsThreadSafe() const override {
        return true;
    }
};

}

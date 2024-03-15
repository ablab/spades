//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"

namespace io {

namespace binary {

template<typename Index>
class BaseCoverageIO : public IOSingle<Index> {
public:
    BaseCoverageIO(const char *name, const char *ext):
            IOSingle<Index>(name, ext) {
    }

    void SaveImpl(BinOStream &str, const Index &index) override {
        for (auto e : index.g().canonical_edges()) {
            str << e << index.RawCoverage(e);
        }
        str << uint64_t(0);  // null-term
    }

    void LoadImpl(BinIStream &str, Index &index) override {
        while (true) {
            uint64_t eid;
            str >> eid;
            if (!eid) break; // null-term
            uint32_t cov;
            str >> cov;
            index.SetRawCoverage(eid, cov);
            auto conj_eid = index.g().conjugate(typename Index::EdgeId(eid));
            if (eid != conj_eid.int_id()) {
                index.SetRawCoverage(conj_eid, cov);
            }
        }
    }
};

template<typename Graph>
class CoverageIO : public BaseCoverageIO<omnigraph::CoverageIndex<Graph>> {
public:
    typedef BaseCoverageIO<omnigraph::CoverageIndex<Graph>> base;
    CoverageIO():
            base("coverage", ".cvr") {}
};

template<typename Graph>
struct IOTraits<omnigraph::CoverageIndex<Graph>> {
    typedef CoverageIO<Graph> Type;
};

template<typename Graph>
class FlankingCoverageIO : public BaseCoverageIO<omnigraph::FlankingCoverage<Graph>> {
public:
    typedef BaseCoverageIO<omnigraph::FlankingCoverage<Graph>> base;
    FlankingCoverageIO():
            base("flanking coverage", ".flcvr") {}
};

template<typename Graph>
struct IOTraits<omnigraph::FlankingCoverage<Graph>> {
    typedef FlankingCoverageIO<Graph> Type;
};

} // namespace binary

} //namespace io

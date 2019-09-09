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
        size_t edges_cnt = index.g().e_size();

        str << edges_cnt;
        for (auto it = index.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto e = *it;
            str << e.int_id() << index.RawCoverage(e);
        }
    }

    void LoadImpl(BinIStream &str, Index &index) override {
        size_t edges_cnt;
        str >> edges_cnt;
        for (size_t i = 0; i < edges_cnt; ++i) {
            uint64_t eid;
            uint32_t cov;
            str >> eid >> cov;
            index.SetRawCoverage(eid, cov);
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

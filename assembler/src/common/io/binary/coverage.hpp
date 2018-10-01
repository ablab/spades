//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"
#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"

namespace io {

namespace binary {

template<typename Index>
class BaseCoverageIO : public IOSingle<Index, EdgeMapper<Index>> {
public:
    typedef EdgeMapper<Index> Mapper;

    BaseCoverageIO(const char *name, const char *ext):
            IOSingle<Index, Mapper>(name, ext) {
    }

private:
    void SaveImpl(BinSaveFile &file, const Index &index) override {
        for (auto it = index.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto e = *it;
            file << e.int_id() << index.RawCoverage(e);
        }
    }

    void LoadImpl(BinLoadFile &file, Index &index, const Mapper &mapper) override {
        size_t e;
        while (file >> e) { //Read until the end
            auto eid = mapper[e];
            auto cov = file.Read<unsigned>();
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

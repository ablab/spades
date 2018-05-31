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

template<typename Index>
class BaseCoverageIO : public IOSingle<Index> {
public:
    typedef IdMapper<typename Index::EdgeId> Mapper;
    BaseCoverageIO(const Mapper &mapper, const char *name, const char *ext):
            IOSingle<Index>(name, ext), mapper_(mapper) {}
private:
    void SaveImpl(SaveFile &file, const Index &index) {
        for (auto it = index.g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            auto e = *it;
            file << e.int_id() << index.RawCoverage(e);
        }
    }

    void LoadImpl(LoadFile &file, Index &index) {
        size_t e;
        while (file >> e) { //Read until the end
            auto eid = mapper_[e];
            auto cov = file.Read<unsigned>();
            index.SetRawCoverage(eid, cov);
        }
    }

    const Mapper &mapper_;
};

template<typename Graph>
class CoverageIO : public BaseCoverageIO<omnigraph::CoverageIndex<Graph>> {
public:
    typedef BaseCoverageIO<omnigraph::CoverageIndex<Graph>> base;
    CoverageIO(const typename base::Mapper &mapper):
            base(mapper, "coverage", ".cvr") {}
};

template<typename Graph>
class FlankingCoverageIO : public BaseCoverageIO<omnigraph::FlankingCoverage<Graph>> {
public:
    typedef BaseCoverageIO<omnigraph::FlankingCoverage<Graph>> base;
    FlankingCoverageIO(const typename base::Mapper &mapper):
            base(mapper, "flanking coverage", ".flcvr") {}
};

}

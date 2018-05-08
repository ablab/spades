//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "io/id_mapper.hpp"

namespace io {

template<typename Index>
class CoverageIO : public IOSingle<Index> {
public:
    typedef IdMapper<typename Index::EdgeId> Mapper;
    CoverageIO(const Mapper &mapper):
            IOSingle<Index>("coverage", ".cvr"), mapper_(mapper) {}
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

}

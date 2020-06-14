//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/adt/pack.hpp"

#include <boost/noncopyable.hpp>

namespace debruijn_graph {

class GraphPack: public adt::pack, private boost::noncopyable {
public:
    GraphPack(size_t k_,
               const std::string &workdir_, size_t lib_count,
               const std::vector<std::string> &genome = std::vector<std::string>(0),
               size_t flanking_range = 50,
               size_t max_mapping_gap = 0,
               size_t max_gap_diff = 0,
               bool detach_indices = true);

    size_t k() const { return k_; }
    const std::string &workdir() const { return workdir_; }

    void FillQuality();
    void ClearQuality();

    void EnsureIndex();
    void EnsureBasicMapping();
    void EnsureQuality();
    void EnsurePos();
    void EnsureDebugInfo();

    void InitRRIndices();
    void ClearRRIndices();
    void ClearPaths();

    void DetachAll();
    void DetachEdgeIndex();

    void PrepareForStage(const char*);

private:
    size_t k_;
    std::string workdir_;
};

} // namespace debruijn_graph

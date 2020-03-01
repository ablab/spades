//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "long_read_mapper.hpp"

namespace debruijn_graph {

PathWithMappingInfo::PathWithMappingInfo(std::vector<EdgeId> && path, Range && range) 
    : Path_(std::move(path))
    , MappingRangeOntoRead_(range)
{}

std::mutex PathsWithMappingInfoStorageStorageLock;
std::vector<std::unique_ptr<path_extend::BidirectionalPath>> PathsWithMappingInfoStorageStorage;

} // namespace debruijn_graph


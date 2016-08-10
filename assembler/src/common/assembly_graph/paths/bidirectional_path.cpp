//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * bidirectional_path.cpp
 *
 *  Created on: Jun 25, 2015
 *      Author: andrey
 */

#include "utils/standard_base.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"

namespace path_extend {

std::atomic<uint64_t> BidirectionalPath::path_id_{0};

}

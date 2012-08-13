//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * standart.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard_base.hpp"

// utils
#include "cpp_utils.hpp"
#include "path_helper.hpp"

#include "simple_tools.hpp"

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"

// omni
#include "omni/paired_info.hpp"
#include "omni/total_labeler.hpp"

#include "runtime_k.hpp"
// common typedefs

namespace debruijn_graph
{
    // todo: rename
    typedef io::IReader<io::SingleRead> SingleReadStream;
    typedef io::IReader<io::PairedRead> PairedReadStream;
    typedef io::ConvertingReaderWrapper UnitedStream;
} // namespace debruijn_graph


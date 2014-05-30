//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard_base.hpp"

// log
#include "logger/logger.hpp"

// utils
#include "cpp_utils.hpp"
#include "path_helper.hpp"

#include "simple_tools.hpp"

// longseq
#include "longseq.hpp"

// config
#include "cap_config_struct.hpp"

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include "io/multifile_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/osequencestream.hpp"

namespace cap {
typedef io::SingleRead Contig;
typedef io::ReadStream<Contig> ContigStream;
typedef std::shared_ptr<ContigStream> ContigStreamPtr;
typedef	io::MultifileStream<io::SingleRead> CompositeContigStream;
typedef	io::RCWrapper<io::SingleRead> RCWrapper;
typedef io::ReadStreamList<Contig> ContigStreams;
}

// debruijn
#include "debruijn_graph.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"

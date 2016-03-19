//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "dev_support/standard_base.hpp"

// log
#include "logger/logger.hpp"

// utils
#include "dev_support/cpp_utils.hpp"
#include "io_impl/path_helper.hpp"

#include "dev_support/simple_tools.hpp"

// longseq
#include "longseq.hpp"

// config
#include "cap_config_struct.hpp"

// io
#include "io/reads_io/ireader.hpp"
#include "io/reads_io/converting_reader_wrapper.hpp"
#include "io/reads_io/vector_reader.hpp"
#include "io/reads_io/multifile_reader.hpp"
#include "io/reads_io/rc_reader_wrapper.hpp"
#include "io/reads_io/osequencestream.hpp"

namespace cap {
typedef io::SingleRead Contig;
typedef io::ReadStream<Contig> ContigStream;
typedef std::shared_ptr<ContigStream> ContigStreamPtr;
typedef	io::MultifileStream<io::SingleRead> CompositeContigStream;
typedef	io::RCWrapper<io::SingleRead> RCWrapper;
typedef io::ReadStreamList<Contig> ContigStreams;
}

// debruijn
#include "assembly_graph/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "graph_construction.hpp"

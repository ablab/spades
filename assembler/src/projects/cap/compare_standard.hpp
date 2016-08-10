//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"

// log
#include "utils/logger/logger.hpp"

// utils
#include "utils/cpp_utils.hpp"
#include "utils/path_helper.hpp"

#include "utils/simple_tools.hpp"

// longseq
#include "longseq.hpp"

// config
#include "cap_config_struct.hpp"

// io
#include "io/reads/ireader.hpp"
#include "io/reads/converting_reader_wrapper.hpp"
#include "io/reads/vector_reader.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "io/reads/osequencestream.hpp"

namespace cap {
typedef io::SingleRead Contig;
typedef io::ReadStream<Contig> ContigStream;
typedef std::shared_ptr<ContigStream> ContigStreamPtr;
typedef    io::MultifileStream<io::SingleRead> CompositeContigStream;
typedef    io::RCWrapper<io::SingleRead> RCWrapper;
typedef io::ReadStreamList<Contig> ContigStreams;
}

// debruijn
#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "modules/graph_construction.hpp"

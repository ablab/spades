//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include "io/multifile_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/ofastastream.hpp"

namespace cap {
typedef io::SingleRead Contig;
typedef io::IReader<Contig> ContigStream;
typedef	io::MultifileReader<io::SingleRead> CompositeContigStream;
typedef	io::RCReaderWrapper<io::SingleRead> RCWrapper;
typedef LongSeq<MultiPolynomialHash<3, uint64_t> > LSeq;
}

// debruijn
#include "new_debruijn.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"

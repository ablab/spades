#pragma once

#include "standard_base.hpp"

// log
#include "logger/logger.hpp"

// utils
#include "cpp_utils.hpp"
#include "fs_path_utils.hpp"

#include "simple_tools.hpp"

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include "io/multifile_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/ofastastream.hpp"

namespace compare {
typedef io::SingleRead Contig;
typedef io::IReader<Contig> ContigStream;
typedef	io::MultifileReader<io::SingleRead> CompositeContigStream;
typedef	io::RCReaderWrapper<io::SingleRead> RCWrapper;
}

// debruijn
#include "new_debruijn.hpp"
#include "graph_pack.hpp"
#include "graph_construction.hpp"

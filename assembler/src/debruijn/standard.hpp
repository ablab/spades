/*
 * standart.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once
#include "standard_base.hpp"
#include "k.hpp"

//==our
// log
#include "logger/logger.hpp"

// utils
#include "cpp_utils.hpp"
#include "fs_path_utils.hpp"

#include "simple_tools.hpp"

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"

// omni
#include "omni/paired_info.hpp"
#include "omni/total_labeler.hpp"

// common typedefs

namespace debruijn_graph
{
    // todo: rename
    typedef io::IReader<io::SingleRead> SingleReadStream;
    typedef io::IReader<io::PairedRead> PairedReadStream;
    typedef io::ConvertingReaderWrapper UnitedStream;
} // namespace debruijn_graph

inline bool make_dir(fs::path p)
{
	namespace fs = boost::filesystem;
	if (fs::is_directory(p) || fs::create_directories(p))
		return true;

	WARN("Can't create directory " << p);
	return false;
}

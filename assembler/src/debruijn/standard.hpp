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
// utils
#include "cpp_utils.hpp"
#include "fs_path_utils.hpp"

// err handling
#include "stacktrace.hpp"

#ifndef NDEBUG
namespace boost
{
inline void assertion_failed(char const * expr, char const * function, char const * file, long line)
{
        std::cerr << "Aborted by assert: " << std::endl;
        print_stacktrace();
        __assert_fail (expr, file, line, function);
}
} // namespace boost
#endif


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

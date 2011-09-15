/*
 * standart.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

//==crt and stl
#include <cstdlib>
#include <time.h>
#include <signal.h>
#include <execinfo.h>

#include <algorithm>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <utility>

using std::list;
using std::map;
using std::vector;
using std::set;
using std::make_pair;

//==boost
#include <boost/algorithm/string.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>

#include <boost/filesystem.hpp>

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

using boost::bind;
using boost::function;

using boost::shared_ptr;
using boost::scoped_ptr;
using boost::make_shared;

namespace fs = boost::filesystem;

using boost::optional;
using boost::in_place;

//==sys
#include <sys/stat.h>
#include <sys/types.h>


//==our
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

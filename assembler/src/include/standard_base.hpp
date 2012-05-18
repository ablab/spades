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

//==crt and stl
#include <memory>
#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <signal.h>
#include <execinfo.h>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <array>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <stack>
#include <cmath>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::map;
using std::multimap;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::array;
using std::set;
using std::string;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::min;
using std::max;
using std::abs;
using std::stringstream;
using std::numeric_limits;
using std::deque;
using std::ostream_iterator;
using std::stack;
using std::copy;

//==boost

#ifndef NDEBUG
#define BOOST_ENABLE_ASSERT_HANDLER
#endif

#include <boost/assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/bimap.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/smart_ptr.hpp>
#include <boost/make_shared.hpp>

#include <boost/filesystem.hpp>

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

#include <boost/format.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/noncopyable.hpp>
#include <boost/algorithm/string.hpp>

using boost::bind;
using boost::ref;
using boost::cref;
using boost::function;

using boost::bimap;
using boost::shared_ptr;
using boost::scoped_ptr;
using boost::make_shared;

namespace fs = boost::filesystem;

using boost::optional;
using boost::none;
using boost::in_place;

using boost::format;

using boost::lexical_cast;
using boost::noncopyable;

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

//==sys
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

//our
#include "xmath.h"

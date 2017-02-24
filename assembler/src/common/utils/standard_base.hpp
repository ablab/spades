//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
#include <map>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <cmath>
#include <limits>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
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
using std::ostream_iterator;
using std::copy;

using std::shared_ptr;
using std::make_shared;

//==boost

#ifndef NDEBUG
#define BOOST_ENABLE_ASSERT_HANDLER
#endif

#include <boost/optional.hpp>

#include <boost/noncopyable.hpp>

using boost::optional;
using boost::make_optional;
using boost::none;

using boost::noncopyable;

// err handling
#include "utils/stacktrace.hpp"

// path manipulation instead of boost filesystem
#include "utils/path_helper.hpp"
using path::make_dir;
using path::remove_dir;

#ifndef NDEBUG
namespace boost {
inline void assertion_failed(char const * expr, char const * function,
                             char const * file, long line) {
  std::cerr << "Aborted by assert: " << std::endl;
  print_stacktrace();
#if __DARWIN_UNIX03
  __assert_rtn (expr, file, (int)line, function);
#elif __DARWIN
  __assert (expr, file, (int)line, function);
#else
  __assert_fail (expr, file, (unsigned)line, function);
#endif
}

inline void assertion_failed_msg(char const * expr, char const * msg,
                                 char const * function, char const * file,
                                 long line) {
  std::cerr << "Aborted by assert: " << msg << std::endl;
  print_stacktrace();
#if __DARWIN_UNIX03
  __assert_rtn (expr, file, (int)line, function);
#elif __DARWIN
  __assert (expr, file, (int)line, function);
#else
  __assert_fail (expr, file, (unsigned)line, function);
#endif
}

} // namespace boost

#endif // NDEBUG

//==sys
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

//our
//math
#include "math/xmath.h"
#include "func/func.hpp"
#include "utils/verify.hpp"
// log
#include "utils/logger/logger.hpp"



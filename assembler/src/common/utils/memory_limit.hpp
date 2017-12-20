//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdlib>

namespace utils {

void limit_memory(size_t limit);
size_t get_memory_limit();
size_t get_max_rss();
size_t get_used_memory();
size_t get_free_memory();

}

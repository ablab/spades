/*
 * include.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: lab42
 */
#pragma once

// WTF: This header should never exist. Include only WHAT YOURE USING, NOT THE WHOLE PROJECT.

#include "config_struct.hpp"
#include "standard.hpp"
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"
#include "utils.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "io/osequencestream.hpp"
#include "io/ireader.hpp"
#include <vector>

#include <string>
// WTF: Use enum
#define MAX_VARIANTS 7
#define DELETION 5
#define INSERTION 6
#define UNDEFINED 4
namespace corrector {
const char pos_to_var[MAX_VARIANTS]= {'A','C','G','T','N','D','I'};
const set<char> valid_variant = {'a', 'A', 'c', 'C', 'g', 'G', 't', 'T', 'n', 'N', 'I', 'D'};
const int var_to_pos[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


// WTF: This is very inefficient. Make static table
inline static bool IsValidVariant(char C) {
	return valid_variant.find(C) != valid_variant.end();
}

};

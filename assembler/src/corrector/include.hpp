/*
 * include.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: lab42
 */
#pragma once
#include "standard.hpp"
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "io/osequencestream.hpp"
#include <string>
#define MAX_VARIANTS 7
#define DELETION 5
#define INSERTION 6
#define UNDEFINED 4


//const map<char, int> nt_to_pos = {{'a', 0}, {'A', 0}, {'c', 1}, {'C', 1}, {'g', 2}, {'G', 2}, {'t', 3}, {'T', 3}, {'D', 4}, {'I', 5}};
//const map< char, char> pos_to_nt = {{0, 'A'},  {1, 'C'},  {2, 'G'}, {3, 'T'}, {4, 'N'} {5, 'D'}, {6, 'I'}};
const char pos_to_var[MAX_VARIANTS]= {'A','C','G','T','N','D','I'};
const set<char> valid_variant = {'a', 'A', 'c', 'C', 'g', 'G', 't', 'T', 'n', 'N', 'I', 'D'};
const int var_to_pos[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
inline static bool IsValidVariant(char C) {
	return valid_variant.find(C) != valid_variant.end();
}

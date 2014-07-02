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
#define MAX_VOTES 6
#define DELETION 4
#define INSERTION 5


const map<char, int> nt_to_pos = {{'a', 0}, {'A', 0}, {'c', 1}, {'C', 1}, {'g', 2}, {'G', 2}, {'T', 3}, {'T', 3}, {'D', 4}, {'I', 5}};
const map< char, char> pos_to_nt = {{0, 'A'},  {1, 'C'},  {2, 'G'}, {3, 'T'}, {4, 'D'}, {5, 'I'}};

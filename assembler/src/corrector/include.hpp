/*
 * include.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: lab42
 */
#pragma once

// WTF: Why it's here? Include *only* what you're using
#include "standard_base.hpp"

// WTF: Make it out of sync: introduce the last state (with a value) to determine the max amount
#define MAX_VARIANTS 7
enum Variants {NuclA, NuclC, NuclT, NuclG, Undefined, Deletion, Insertion};

namespace corrector {
// WTF: constexpr
const char pos_to_var[MAX_VARIANTS] = { 'A', 'C', 'G', 'T', 'N', 'D', 'I' };
const size_t var_to_pos[128] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
}
;

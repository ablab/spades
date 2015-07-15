//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * include.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: lab42
 */
#pragma once

#include <sys/types.h>

enum Variants {
    NuclA = 0,
    NuclC = 1,
    NuclT = 2,
    NuclG = 3,
    Undefined = 4,
    Deletion = 5,
    Insertion = 6,
    VariantsNumber = 7
};
#define MAX_VARIANTS Variants::VariantsNumber

namespace corrector {
constexpr char pos_to_var[MAX_VARIANTS] = { 'A', 'C', 'G', 'T', 'N', 'D', 'I' };
constexpr size_t var_to_pos[128] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 2, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
}
;

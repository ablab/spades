//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <map>
#include <set>
#include <iostream>

#include "math.h"

using namespace std;
using namespace io;

namespace dipspades {

bool is_intersection_exist(Range r1, Range r2){
    if(r1.end_pos <= r2.start_pos || r2.end_pos <= r1.start_pos)
        return false;
    else
        return true;
}

Range get_intersection_of_ranges(Range r1, Range r2){
    VERIFY(is_intersection_exist(r1, r2));
    size_t max_start = max<size_t>(r1.start_pos, r2.start_pos);
    size_t min_end = min<size_t>(r1.end_pos, r2.end_pos);

    Range r(max_start, min_end);
    return r;
}

pair<size_t, size_t> project_init_range_to_new(Range old_map_rg, Range new_map_rg, Range old_init_rg){

    size_t start_pos, end_pos;

    int shift_start = int(new_map_rg.start_pos) - int(old_map_rg.start_pos);
    double start_coeff = double(shift_start) / double(old_map_rg.end_pos - old_map_rg.start_pos);
    start_pos = old_init_rg.start_pos + int(start_coeff * double(old_init_rg.end_pos - old_init_rg.start_pos));

    int shift_end = int(new_map_rg.end_pos) - int(old_map_rg.end_pos);
    double end_coeff = double(shift_end) / double(old_map_rg.end_pos - old_map_rg.start_pos);
    end_pos = old_init_rg.end_pos + int(end_coeff * double(old_init_rg.end_pos - old_init_rg.start_pos));

    return pair<size_t, size_t>(start_pos, end_pos);
}

bool is_range_pair_correct(pair<size_t, size_t> p){
    return p.first < p.second;
}

}

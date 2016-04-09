//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "positional_read.hpp"

#include <sstream>
using namespace std;

namespace corrector {

string position_description::str() const {
    stringstream ss;
    for (int i = 0; i < MAX_VARIANTS; i++) {
        ss << pos_to_var[i];
        ss << ": " << votes[i] << "; ";
    }
    return ss.str();
}

void position_description::clear() {
    for (size_t i = 0; i < MAX_VARIANTS; i++) {
        votes[i] = 0;
        insertions.clear();
    }
}
};

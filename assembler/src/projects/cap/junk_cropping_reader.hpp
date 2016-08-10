//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "utils/standard_base.hpp"
#include "io/reads/delegating_reader_wrapper.hpp"

namespace cap {

class JunkCroppingWrapper : public io::DelegatingWrapper<io::SingleRead> {
    typedef io::DelegatingWrapper<io::SingleRead> base;
    vector<pair<size_t, size_t>> coordinates_ladder_;

    bool IsGoodSymbol(char c) {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    }

public:
    JunkCroppingWrapper(base::ReadStreamPtrT reader) : base(reader) {

    }

    JunkCroppingWrapper& operator>>(io::SingleRead& read) {
        base::operator >>(read);
        coordinates_ladder_.clear();
        string orig_string = read.GetSequenceString();
        string orig_qual = read.GetQualityString();
        string cropped = "";
        string cropped_qual = "";
        coordinates_ladder_.push_back(make_pair(0, 0));
        for (size_t coord = 0; coord < orig_string.size(); ++coord) {
            if (coord > 0 && (IsGoodSymbol(orig_string[coord - 1]) ^ IsGoodSymbol(orig_string[coord]))) {
                coordinates_ladder_.push_back(make_pair(cropped.size(), coord));
            }
            if (IsGoodSymbol(orig_string[coord])) {
                cropped += orig_string[coord];
                cropped_qual += orig_qual[coord];
            }
        }
        coordinates_ladder_.push_back(make_pair(cropped.size(), orig_string.size()));
        read = io::SingleRead(read.name(), cropped, cropped_qual);
        return *this;
    }

    vector<pair<size_t, size_t>> coordinates_ladder() const {
        return coordinates_ladder_;
    }
};

}

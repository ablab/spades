//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/library_fwd.hpp"
#include <functional>

namespace io {

static inline std::pair<bool, bool> GetRCFlags(LibraryOrientation orientation) {
    switch (orientation) {
        case LibraryOrientation::RR:
            return {true, true};
        case LibraryOrientation::FR:
            return {false, true};
        case LibraryOrientation::RF:
            return {true, false};
        default:
            return {false, false};
    }
}

class OrientationChanger {
  public:
    OrientationChanger(LibraryOrientation orientation) {
        std::tie(rc_left_, rc_right_) = GetRCFlags(orientation);
    }

    template<class ReadType>
    void operator()(ReadType &read) const {
        // FIXME: do this inplace
        if (rc_left_)
            read.first() = !read.first();
        if (rc_right_)
            read.second() = !read.second();
    }

  private:
    bool rc_left_ = false;
    bool rc_right_ = false;
};

}

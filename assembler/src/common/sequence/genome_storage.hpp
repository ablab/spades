//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include "sequence.hpp"
#include "nucl.hpp"

class GenomeStorage {
    std::string s_;
public:
    GenomeStorage() {
    }

    GenomeStorage(const std::string &s): s_(s) {
    }

    //TODO exterminate this where possible
    Sequence GetSequence() const {
        stringstream ss;
        size_t l = 0, r = 0;
        for(size_t i = 0; i < s_.size(); i++) {
            if (!is_nucl(s_[i]) ) {
                if (r > l) {
                    ss << s_.substr(l, r - l);
                }
                r = i + 1;
                l = i + 1;
            } else {
                r++;
            }
        }
        if (r > l) {
            ss << s_.substr(l, r - l);
        }
        return Sequence(ss.str());
    }

    void SetSequence(const Sequence &s) {
        s_ = s.str();
    }

    std::string str() const {
        return s_;
    }

    size_t size() const {
        return s_.size();
    }
};


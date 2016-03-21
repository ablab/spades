//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by lab42 on 8/19/15.
//

#include "genome_storage.hpp"
#include "data_structures/sequence/nucl.hpp"
using namespace std;

namespace debruijn_graph {
//TODO exterminate this where possible
    Sequence GenomeStorage::GetSequence() const{
        stringstream ss;
        size_t l = 0, r = 0;
        for(size_t i = 0; i < s_.size(); i++) {
            if (! is_nucl(s_[i]) ) {
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
    void GenomeStorage::SetSequence(const Sequence &s) {
        s_ = s.str();
    }
    string GenomeStorage::str() const{
        return s_;
    }
    size_t GenomeStorage::size() const {
        return s_.size();
    }
}
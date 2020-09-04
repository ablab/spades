//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once


namespace kmer {

template<class Seq>
struct KMerBucketPolicy {
    typedef typename Seq::hash hash;

    unsigned operator()(const Seq &s, unsigned total) const {
        return (unsigned)(hash()(s) % total);
    }
};
    

}

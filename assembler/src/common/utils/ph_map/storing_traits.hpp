//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"

namespace utils {

struct SimpleStoring {
    template<class K, class Container>
    static auto get_value(const Container &values, const K& key) {
        return values[key.idx()];
    }

    template<class K, class V, class Container>
    static void set_value(Container &values, const K& key, const V& value) {
        values[key.idx()] = value;
    }

    static constexpr bool IsInvertable() {
        return false;
    }
};

struct InvertableStoring {
    struct default_inverter {
        template<typename K, typename V>
        V operator()(const V& v, const K& k) const {
            return v.conjugate(k);
        }
    };

    struct trivial_inverter {
        template<typename K, typename V>
        V operator()(const V& v, const K& /*k*/) const {
            return v;
        }
    };

    template<class K, class Container, class F = default_inverter>
    static auto get_value(const Container &values, const K& key,
                          const F& inverter = F()) {
        if (key.is_minimal())
            return values[key.idx()];
        else
            return inverter(values[key.idx()], key);
    }

    template<class K, class V, class Container, class F = default_inverter>
    static void set_value(Container& values, const K& key, const V& value,
                          const F& inverter = F()) {
        VERIFY(key.idx() < values.size());
        if (key.is_minimal()) {
            values[key.idx()] = value;
        } else {
            values[key.idx()] = inverter(value, key);
        }
    }

    static constexpr bool IsInvertable() {
        return true;
    }
};

typedef InvertableStoring DefaultStoring;

//FIXME make functor
template<class StoringType>
struct StoringTypeFilter {
};

template<>
struct StoringTypeFilter<SimpleStoring> {
    template<class Kmer>
    bool operator()(const Kmer &/*kmer*/) const {
        return true;
    }

    template<class Kmer>
    static bool filter(const Kmer &/*kmer*/) {
        return true;
    }
};

template<>
struct StoringTypeFilter<InvertableStoring> {
    template<class Kmer>
    bool operator()(const Kmer &kmer) const {
        return filter(kmer);
    }

    template<class Kmer>
    static bool filter(const Kmer &kmer) {
        return kmer.IsMinimal();
    }
};

}

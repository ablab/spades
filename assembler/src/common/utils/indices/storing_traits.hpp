//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
/*
 * key_with_hash.hpp
 *
 *  Created on: Nov 7, 2013
 *      Author: anton
 */

#include "values.hpp"

namespace debruijn_graph {


struct SimpleStoring {
    template<class K, class V>
    static V get_value(const ValueArray<V> &values, const K& key) {
        return values[key.idx()];
    }

    template<class K, class V>
    static void set_value(ValueArray<V> &values, const K& key, const V& value) {
        values[key.idx()] = value;
    }

    static bool IsInvertable() {
        return false;
    }
};

struct InvertableStoring {
    template<typename V>
    struct default_inverter {
        template<typename K>
        V operator()(const V& v, const K& k) const {
            return v.conjugate(k);
        }
    };

    template<typename V>
    struct trivial_inverter {
        template<typename K>
        V operator()(const V& v, const K& /*k*/) const {
            return v;
        }
    };

    template<class K, class V, class F = default_inverter<V>>
    static V get_value(const ValueArray<V> &values, const K& key,
                       const F& inverter = F()) {
        if (key.is_minimal())
            return values[key.idx()];
        else
            return inverter(values[key.idx()], key);
    }

    template<class K, class V, class F = default_inverter<V>>
    static void set_value(ValueArray<V>& values, const K& key, const V& value,
                          const F& inverter = F()) {
        VERIFY(key.idx() < values.size());
        if (key.is_minimal()) {
            values[key.idx()] = value;
        } else {
            values[key.idx()] = inverter(value, key);
        }
    }

    static bool IsInvertable() {
        return true;
    }
};

typedef InvertableStoring DefaultStoring;

}

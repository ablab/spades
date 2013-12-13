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
    static void  set_value(ValueArray<V> &values, const K& key, const V& value) {
        values[key.idx()] = value;
    }

    static bool IsInvertable() {
        return false;
    }
};

struct InvertableStoring {
    template<class K, class V>
    static V get_value(const ValueArray<V> &values, const K& key) {
        if(key.is_minimal())
            return values[key.idx()];
        else
            return values[key.idx()].conjugate(key);
    }

    template<class K, class V>
    static void set_value(ValueArray<V> &values, const K& key, const V& value) {
        if(key.is_minimal())
            values[key.idx()] = value;
        else
            values[key.idx()] = value.conjugate(key);
    }

    static bool IsInvertable() {
        return true;
    }
};

typedef InvertableStoring DefaultStoring;

}

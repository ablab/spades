//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace omnigraph {

template<class T>
class PairedElementManipulationHelper {
public:
    bool IsMinimal(T t) const {
        return !(t->conjugate_ < t);
    }

    T MinimalFromPair(T t) const {
        if (IsMinimal(t)) {
            return t;
        } else {
            return t->conjugate_;
        }
    }

    T& GetElementToManipulate(T t) const {
        return t->conjugate_;
    }

    T& ToManipulateFromPair(T t) const {
        return GetElementToManipulate(MinimalFromPair(t));
    }
};

template<class T>
class GraphElementLock : PairedElementManipulationHelper<T> {
    PairedElementManipulationHelper<T> helper_;
    restricted::PurePtrLock<T> inner_lock_;

public:
    GraphElementLock(T  t) :
        inner_lock_(helper_.ToManipulateFromPair(t))
    {
    }

};

/**
 * Do not use with locks on same graph elements!
 */
template<class T>
class GraphElementMarker {
    PairedElementManipulationHelper<T> helper_;
    restricted::PurePtrMarker<T> marker_;
public:

    void mark(T t) {
        marker_.mark(helper_.ToManipulateFromPair(t));
    }

    void unmark(T t) {
        marker_.unmark(helper_.ToManipulateFromPair(t));
    }

    bool is_marked(T t) const {
        return marker_.is_marked(helper_.ToManipulateFromPair(t));
    }
};
}

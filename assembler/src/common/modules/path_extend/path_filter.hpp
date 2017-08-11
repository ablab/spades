//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * path_filter.hpp
 *
 *  Created on: Mar 14, 2012
 *      Author: andrey
 */

#ifndef PATH_FILTER_HPP_
#define PATH_FILTER_HPP_

#include "assembly_graph/paths/bidirectional_path.hpp"

namespace path_extend {

typedef func::AbstractPredicate<const BidirectionalPath&> AbstractPathCondition;

class EmptyPathCondition: public AbstractPathCondition {
public:
    EmptyPathCondition() {}

    bool Check(checked_type p) const override {
        return p.Empty();
    }
};

class LengthPathCondition: public AbstractPathCondition {
    size_t min_length_;
public:
    LengthPathCondition(size_t min_length): min_length_(min_length) {}

    bool Check(checked_type p) const override {
        return p.Length() <= min_length_;
    }
};

class CoveragePathCondition: public AbstractPathCondition {
    const Graph& g_;
    double cov_;

public:
    CoveragePathCondition(const Graph& g, double cov): g_(g), cov_(cov) {}

    bool Check(checked_type p) const override {
        for (size_t i = 0; i < p.Size(); ++i) {
            if (math::gr(g_.coverage(p[i]), cov_))
                return false;
        }
        return true;
    }
};

class IsolatedPathCondition: public AbstractPathCondition {
    const Graph& g_;
public:
    IsolatedPathCondition(const Graph& g): g_(g) {}

    bool Check(checked_type p) const override {
        if (p.Empty())
            return true;

        if (p.Size() <= 2) {
            auto v1 = g_.EdgeStart(p.Front());
            auto v2 = g_.EdgeEnd(p.Back());

            return g_.IncomingEdgeCount(v1) == 0 &&
                   g_.OutgoingEdgeCount(v2) == 0;
        }
        return false;
    }
};

}

#endif /* PATH_FILTER_HPP_ */

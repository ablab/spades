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

class CopyOnWritePathFilter {

protected:
    const Graph& g;

public:
    CopyOnWritePathFilter(const Graph& g_): g(g_) {
    }

    virtual bool predicate(BidirectionalPath& path) = 0;

    PathContainer filter(PathContainer& paths) {
        PathContainer result;

        for (size_t i = 0; i < paths.size(); ++i) {
            if (predicate(*paths.Get(i)) || predicate(*paths.GetConjugate(i))) {
                result.AddPair(paths.Get(i), paths.GetConjugate(i));
            }
        }

        return result;
    }

};


class IdFilter: public CopyOnWritePathFilter {

protected:
    std::set<size_t> ids;

public:

    IdFilter(const Graph& g_, std::set<size_t> ids_): CopyOnWritePathFilter(g_), ids(ids_) {
    }

    virtual bool predicate(BidirectionalPath& path) {
        return ids.count(path.GetId()) > 0;
    }
};


class DuplicateFilter {

protected:
    const Graph& g;

public:
    DuplicateFilter(const Graph& g_): g(g_) {
    }

    PathContainer filter(PathContainer& paths) {
        PathContainer result;

        for (size_t i = 0; i < paths.size(); ++i) {
            bool duplicate = false;
            for (size_t j = 0; j < result.size(); ++j) {
                if (result[j] == paths[j])
                    duplicate = true;
            }
            if (!duplicate) {
                result.AddPair(paths.Get(i), paths.GetConjugate(i));
            }
        }

        return result;
    }

};

class ErasingPathFilter {

protected:
    const Graph& g;

public:
    ErasingPathFilter(const Graph& g_): g(g_) {
    }

    virtual bool predicate(BidirectionalPath& path) = 0;

    void filter(PathContainer& paths) {
        for (PathContainer::Iterator iter = paths.begin(); iter != paths.end(); ) {
            if (predicate(*iter.get()) || predicate(*iter.getConjugate())) {
                iter = paths.erase(iter);
            }
            else {
                ++iter;
            }
        }
    }

};


class CoveragePathFilter: public ErasingPathFilter {

protected:
    double minCoverage;

public:
    CoveragePathFilter(Graph& g_, double cov): ErasingPathFilter(g_), minCoverage(cov) {

    }

    virtual bool predicate(BidirectionalPath& path) {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (math::ls(g.coverage(path[i]), minCoverage)) {
                return true;
            }
        }
        return false;
    }
};


class LengthPathFilter: public ErasingPathFilter {

protected:
    size_t minLength;

public:
    LengthPathFilter(const Graph& g_, size_t len): ErasingPathFilter(g_), minLength(len) {
    }

    virtual bool predicate(BidirectionalPath& path) {
        return path.Length() <= minLength;
    }
};


class IsolatedPathFilter: public ErasingPathFilter {

protected:
    size_t min_length_;

    double min_cov_;

public:
    IsolatedPathFilter(const Graph& g_, size_t min_length, double min_cov = 10000000.0):
        ErasingPathFilter(g_),
        min_length_(min_length),
        min_cov_(min_cov) {
    }

    virtual bool predicate(BidirectionalPath& path) {
        if (path.Empty())
            return true;

        if (path.Size() <= 2) {
            auto v1 = g.EdgeStart(path.Front());
            auto v2 = g.EdgeEnd(path.Back());

            return g.IncomingEdgeCount(v1) == 0 &&
                g.OutgoingEdgeCount(v2) == 0 &&
                path.Length() < min_length_ &&
                math::ls(path.Coverage(), min_cov_);
        }
        return false;
    }
};

}

#endif /* PATH_FILTER_HPP_ */

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_filter.hpp
 *
 *  Created on: Mar 14, 2012
 *      Author: andrey
 */

#ifndef PATH_FILTER_HPP_
#define PATH_FILTER_HPP_

#include "bidirectional_path.hpp"

namespace path_extend {

class CopyOnWritePathFilter {

protected:
    Graph& g;

public:
    CopyOnWritePathFilter(Graph& g_): g(g_) {
    }

    virtual bool predicate(BidirectionalPath& path) = 0;

    virtual bool conjugateOperator(bool p, bool cp) {
        return p || cp;
    }

    PathContainer filter(PathContainer& paths) {
        PathContainer result;

        for (size_t i = 0; i < paths.size(); ++i) {
            if (conjugateOperator(predicate(*paths.Get(i)), predicate(*paths.GetConjugate(i)))) {
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

    IdFilter(Graph& g_, std::set<size_t> ids_): CopyOnWritePathFilter(g_), ids(ids_) {
    }

    virtual bool predicate(BidirectionalPath& path) {
        return ids.count(path.GetId()) > 0;
    }
};


class ErasingPathFilter {

protected:
    const Graph& g;

public:
    ErasingPathFilter(const Graph& g_): g(g_) {
    }

    virtual bool predicate(BidirectionalPath& path) = 0;

    virtual bool conjugateOperator(bool p, bool cp) {
        return p && cp;
    }

    void filter(PathContainer& paths) {
        for (PathContainer::Iterator iter = paths.begin(); iter != paths.end(); ) {
            if (!conjugateOperator(predicate(*iter.get()), predicate(*iter.getConjugate()))) {
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
                return false;
            }
        }
        return true;
    }
};


class LengthPathFilter: public ErasingPathFilter {

protected:
    size_t minLength;

public:
    LengthPathFilter(const Graph& g_, size_t len): ErasingPathFilter(g_), minLength(len) {
    }

    virtual bool predicate(BidirectionalPath& path) {
        return path.Length() > minLength;
    }
};

}

#endif /* PATH_FILTER_HPP_ */

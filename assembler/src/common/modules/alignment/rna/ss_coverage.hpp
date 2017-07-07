//
// Created by andrey on 22.05.17.
//

#pragma once

#include <assembly_graph/core/graph.hpp>

namespace debruijn_graph {


class SSCoverageStorage {
public:
    typedef std::unordered_map<EdgeId, double> InnerMap;

private:
    const Graph& g_;

    InnerMap storage_;

public:
    SSCoverageStorage(const Graph& g): g_(g), storage_() {}

    double GetCoverage(EdgeId e, bool reverse = false) const {
        if (reverse) {
            e = g_.conjugate(e);
        }

        auto it = storage_.find(e);
        if (it == storage_.end())
            return 0.0;
        return it->second;
    }

    void IncreaseKmerCount(EdgeId e, size_t count, bool add_reverse = false) {
        storage_[e] += (double) count;
        if (add_reverse)
            storage_[g_.conjugate(e)] += (double) count;
    }

    void Clear() {
        storage_.clear();
    }

    void RecalculateCoverage() {
        for(auto& it : storage_) {
            it.second = it.second / double(g_.length(it.first));
        }
    }

    InnerMap::const_iterator begin() const {
        return storage_.begin();
    }

    InnerMap::const_iterator end() const {
        return storage_.end();
    }
};


}
//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace superpath_index {

template <typename EdgeId>
size_t find_subpath(const std::vector<EdgeId> &subpath, const std::vector<EdgeId> &path) {
    const size_t npos = size_t(-1);

    if (path.size() < subpath.size()) {
        return npos;
    }

    for (size_t i = 0; i <= path.size() - subpath.size(); ++i) {
        if (std::equal(subpath.cbegin(), subpath.cend(), path.cbegin() + i)) {
            return i;
        }
    }
    return npos;
}

template <typename EdgeId>
std::vector<size_t> find_subpaths(const std::vector<EdgeId> &subpath, const std::vector<EdgeId> &path) {
    std::vector<size_t> result;

    if (path.size() < subpath.size()) {
        return result;
    }

    for (size_t i = 0; i <= path.size() - subpath.size(); ++i) {
        if (std::equal(subpath.cbegin(), subpath.cend(), path.cbegin() + i)) {
            result.push_back(i);
        }
    }
    return result;
}

template <typename EdgeId>
class SuperpathIndex {
public:
    SuperpathIndex(const std::vector<std::vector<EdgeId>> &superpaths) : superpaths_{superpaths} {
        for (size_t i = 0; i < superpaths_.size(); ++i) {
            for (EdgeId edge : superpaths_[i]) {
                edge2paths_[edge].push_back(i);
            }
        }
    }

    const auto &operator[](size_t i) const { return superpaths_[i]; }
    size_t size() const { return superpaths_.size(); }
    bool empty() const { return superpaths_.empty(); }
    auto cbegin() const { return superpaths_.cbegin(); }
    auto cend() const { return superpaths_.cend(); }

    std::vector<std::pair<size_t, size_t>> query(const std::vector<EdgeId> &q) const {
        std::vector<std::pair<size_t, size_t>> result;
        for (size_t i : candidates(q)) {
            for (size_t pos : find_subpaths(q, superpaths_[i])) {
                result.push_back({i, pos});
            }
        }

        return result;
    }

private:
    std::vector<std::vector<EdgeId>> superpaths_;
    std::unordered_map<EdgeId, std::vector<size_t>> edge2paths_;

    std::vector<size_t> candidates(const std::vector<EdgeId> &query) const {
        std::vector<size_t> result;
        VERIFY(!query.empty());
        for (size_t i = 0; i < query.size(); ++i) {
            const auto &ind = indices(query[i]);
            if (i == 0) {
                result = ind;
            } else {
                auto it =
                    std::set_intersection(result.cbegin(), result.cend(), ind.cbegin(), ind.cend(), result.begin());
                result.resize(it - result.begin());
            }
        }
        return result;
    }

    const std::vector<size_t> &indices(EdgeId e) const {
        const static std::vector<size_t> empty;

        auto it = edge2paths_.find(e);
        if (it != edge2paths_.end()) {
            return it->second;
        } else {
            return empty;
        }
    }
};
}  // namespace superpath_index

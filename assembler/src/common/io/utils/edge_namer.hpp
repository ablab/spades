//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/header_naming.hpp"
#include "id_mapper.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <string>
#include <functional>

namespace io {
template<class Graph>
using EdgeNamingF = std::function<std::string(const Graph&, typename Graph::EdgeId)>;

template<class Graph>
EdgeNamingF<Graph> IdNamingF(const std::string &prefix = "") {
    using EdgeId = typename Graph::EdgeId;
    return [=](const Graph &g, EdgeId e) {
        return io::MakeContigId(g.int_id(e), prefix);
    };
}

template<class Graph>
EdgeNamingF<Graph> BasicNamingF(const std::string &prefix = "EDGE") {
    using EdgeId = typename Graph::EdgeId;
    return [=](const Graph &g, EdgeId e) {
        return io::MakeContigId(g.int_id(e), g.length(e) + g.k(), g.coverage(e), prefix);
    };
}

template<class Graph>
std::string ErrorNamingF(const Graph& g, typename Graph::EdgeId e) {
    FATAL_ERROR("Unknown edge: " << g.int_id(e));
    return "";
}

template<class Graph>
EdgeNamingF<Graph> MapNamingF(const io::IdMapper<std::string> &id_mapper,
                              EdgeNamingF<Graph> fallback = ErrorNamingF<Graph>) {
    using EdgeId = typename Graph::EdgeId;
    return [&id_mapper, fallback](const Graph &g, EdgeId e) {
        if (!id_mapper.count(g.int_id(e)))
            return fallback(g, e);
        return id_mapper[g.int_id(e)];
    };
}

template<class Graph>
class CanonicalEdgeHelper {
    const Graph &g_;
    using EdgeId = typename Graph::EdgeId;
    const EdgeNamingF<Graph> naming_f_;
    const std::string pos_orient_;
    const std::string neg_orient_;
public:

    CanonicalEdgeHelper(const Graph &g,
                        EdgeNamingF<Graph> naming_f = IdNamingF<Graph>(),
                        const std::string& pos_orient = "+",
                        const std::string& neg_orient = "-") :
            g_(g), naming_f_(naming_f),
            pos_orient_(pos_orient), neg_orient_(neg_orient) {
    }

    bool IsCanonical(EdgeId e) const {
        return e <= g_.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : g_.conjugate(e);
    }

    std::string GetOrientation(EdgeId e) const {
        return IsCanonical(e) ? pos_orient_ : neg_orient_;
    }

    std::string EdgeOrientationString(EdgeId e,
                                      const std::string &delim = "") const {
        return naming_f_(g_, Canonical(e)) + delim + GetOrientation(e);
    }

    std::string EdgeString(EdgeId e) const {
        VERIFY(IsCanonical(e));
        return naming_f_(g_, e);
    }
};
}


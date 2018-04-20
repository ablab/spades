//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "fastg_writer.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "common/io/reads/osequencestream.hpp"

#include <set>
#include <string>
#include <sstream>

using namespace io;
using namespace debruijn_graph;

static std::string FormHeader(const std::string &id,
                              const std::set<std::string>& next_ids) {
    std::stringstream ss;
    ss << id;
    if (next_ids.size() > 0) {
        auto delim = ":";
        for (const auto &s : next_ids) {
            ss  << delim << s;
            delim = ",";
        }
    }
    ss << ";";
    return ss.str();
}

void FastgWriter::WriteSegmentsAndLinks() {
    io::OFastaReadStream os(fn_);
    for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        std::set<std::string> next;
        for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(e))) {
            next.insert(extended_namer_.EdgeOrientationString(next_e));
        }
        os << io::SingleRead(FormHeader(extended_namer_.EdgeOrientationString(e), next),
                             graph_.EdgeNucls(e).str());
    }
}


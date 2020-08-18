//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "assembly_graph/core/action_handlers.hpp"
#include "assembly_graph/core/graph.hpp"
#include "toolchain/edge_label_helper.hpp"

#include <sstream>
#include <iostream>
#include <vector>
#include <map>

namespace debruijn_graph {

//TODO support pos values other than -1
class PositionStorage : public omnigraph::GraphActionHandler<Graph> {
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;

    //currently value is always -1
    std::multimap<EdgeId, int> poss_;

    void Insert(EdgeId e, int pos = -1) {
        poss_.insert(std::make_pair(e, pos));
    }

public:
    PositionStorage(const Graph &g) :
            omnigraph::GraphActionHandler<Graph>(g, "PositionStorage") {}

    void HandleDelete(EdgeId e) override {
        poss_.erase(e);
    }

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
        size_t cumm_len = 0;
        for (EdgeId e : old_edges) {
            for (int p : utils::get_all(poss_, e)) {
                Insert(new_edge, (p < 0) ? -1 : int(cumm_len + p));
            }
            cumm_len += g().length(e);
        }
    }

    void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) override {
        VERIFY_MSG(false, "No support");
    }

    //currently does not handle positions properly
    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override {
        if (poss_.count(old_edge) == 0)
            return;
        if (old_edge == g().conjugate(old_edge)) {
            Insert(new_edge1);
            Insert(g().conjugate(new_edge1));
            Insert(new_edge2);
        } else {
            Insert(new_edge1);
            Insert(new_edge2);
        }
    }

    void Save(std::ostream &os, const io::EdgeNamingF<Graph> &naming_f = io::IdNamingF<Graph>()) const {
        io::CanonicalEdgeHelper<Graph> canonical_helper(g(), naming_f);
        for (const auto &e_p : poss_) {
            //TODO think what coordinate to ouput
            os << canonical_helper.EdgeOrientationString(e_p.first, "\t")
               << '\t' << e_p.second << '\n';
        }
    }

    void Load(std::istream &is, const io::EdgeLabelHelper<Graph> &label_helper) {
        std::string s;
        while (std::getline(is, s)) {
            std::istringstream ss(s);
            std::string label;
            ss >> label;
            EdgeId e = label_helper.edge(label);
            CHECK_FATAL_ERROR(e != EdgeId(), "Couldn't find edge with int id " << std::stoi(label) << " in the graph");

            std::string orient;
            ss >> orient;
            CHECK_FATAL_ERROR(orient == "+" || orient == "-", "Invalid orientation");

            //Currently ignored
            int pos;
            ss >> pos;

            e = (orient == "+") ? e : g().conjugate(e);
            Insert(e);
        }
    }

private:
    DECL_LOGGER("PositionStorage");
};

}
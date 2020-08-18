//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "profile_storage.hpp"

namespace debruijn_graph {
namespace coverage_profiles {

void EdgeProfileStorage::HandleDelete(EdgeId e) {
    profiles_.erase(e);
}

void EdgeProfileStorage::HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) {
    RawAbundanceVector total(sample_cnt_, 0);
    for (EdgeId e : old_edges) {
        Add(total, utils::get(profiles_, e));
    }
    profiles_[new_edge] = std::move(total);
}

void EdgeProfileStorage::HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
    RawAbundanceVector total(utils::get(profiles_, edge1));
    Add(total, utils::get(profiles_, edge2));
    profiles_[new_edge] = std::move(total);
}

void EdgeProfileStorage::HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) {
    AbundanceVector abund = profile(old_edge);
    if (old_edge == g().conjugate(old_edge)) {
        RawAbundanceVector raw1 = MultiplyEscapeZero(abund, g().length(new_edge1));
        profiles_[new_edge1] = raw1;
        profiles_[g().conjugate(new_edge1)] = std::move(raw1);
        profiles_[new_edge2] = MultiplyEscapeZero(abund, g().length(new_edge2));
    } else {
        profiles_[new_edge1] = MultiplyEscapeZero(abund, g().length(new_edge1));
        profiles_[new_edge2] = MultiplyEscapeZero(abund, g().length(new_edge2));
    }
}

void EdgeProfileStorage::Save(std::ostream &os, const io::EdgeNamingF<Graph> &edge_namer) const {
    for (auto it = g().ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        os << edge_namer(g(), e) << '\t';
        auto prof = profile(e);
        std::copy(prof.begin(), prof.end(), std::ostream_iterator<double>(os, "\t"));
        os << '\n';
    }
}

void EdgeProfileStorage::Load(std::istream &is,
                              const io::EdgeLabelHelper<Graph> &label_helper,
                              bool check_consistency) {
    std::string s;
    while (std::getline(is, s)) {
        std::istringstream ss(s);
        std::string label;
        ss >> label;
        EdgeId e = label_helper.edge(label);
        auto p = MultiplyEscapeZero(LoadAbundanceVector(ss), g().length(e));
        profiles_[e] = p;
        profiles_[g().conjugate(e)] = std::move(p);
    }

    if (check_consistency) {
        for (auto it = g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            CHECK_FATAL_ERROR(profiles_.count(e) > 0, "Failed to load profile for one of the edges");
        }
    }
}

}
}

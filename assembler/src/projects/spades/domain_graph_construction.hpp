//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/pipeline/stage.hpp"

namespace debruijn_graph {

//todo rename
class DomainGraphConstruction : public spades::AssemblyStage {
  public:
    DomainGraphConstruction()
        : AssemblyStage("Domain Graph Construction", "domain_graph_construction") {}

    void run(conj_graph_pack &gp, const char*);
};

using ContigAlnInfo = std::unordered_map<std::string, std::string>;
class DomainMatcher {
public:
    void MatchDomains(conj_graph_pack &gp, std::vector<std::string> &domain_filenames);
};

}

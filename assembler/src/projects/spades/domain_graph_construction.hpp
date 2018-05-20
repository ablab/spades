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

struct AlnInfo {
    std::string type;
    std::string name;
    unsigned posl, posr;
    std::string seq;
};

using ContigAlnInfo = std::vector<AlnInfo>;

class DomainMatcher {
public:
    ContigAlnInfo MatchDomains(conj_graph_pack &gp);
};

}

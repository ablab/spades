//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "extract_domains.hpp"
#include "domain_matcher.hpp"

namespace debruijn_graph {

void ExtractDomains::run(graph_pack::GraphPack &gp, const char*) {
    nrps::DomainMatcher().MatchDomains(gp, cfg::get().hm->hmm_set, cfg::get().output_dir);
}

}

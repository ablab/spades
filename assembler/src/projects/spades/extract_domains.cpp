//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "extract_domains.hpp"
#include "domain_matcher.hpp"

namespace debruijn_graph {

void ExtractDomains::run(GraphPack &gp, const char*) {
    nrps::DomainMatcher().MatchDomains(gp, cfg::get().hm->hmm_set, cfg::get().output_dir);
}

}

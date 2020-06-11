//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <pipeline/graph_pack.hpp>

#include <string>
#include <vector>

namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
}

namespace nrps {

struct AlnInfo {
    std::string name;
    std::string type;
    std::string desc;
    unsigned posl, posr;
    std::string seq;
};

using ContigAlnInfo = std::vector<AlnInfo>;

class DomainMatcher {
public:
    ContigAlnInfo MatchDomains(debruijn_graph::GraphPack &gp,
                               const std::string &hmm_set,
                               const std::string &output_dir);
};

}

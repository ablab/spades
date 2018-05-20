//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <string>
#include <vector>

namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
template<class Graph> struct graph_pack;
typedef graph_pack<ConjugateDeBruijnGraph> conj_graph_pack;
}

namespace nrps {

struct AlnInfo {
    std::string type;
    std::string name;
    unsigned posl, posr;
    std::string seq;
};

using ContigAlnInfo = std::vector<AlnInfo>;

class DomainMatcher {
public:
    ContigAlnInfo MatchDomains(debruijn_graph::conj_graph_pack &gp,
                               const std::string &hmm_set,
                               const std::string &output_dir);
};

}

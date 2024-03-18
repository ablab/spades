//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <pipeline/graph_pack.hpp>

#include <filesystem>
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
    ContigAlnInfo MatchDomains(graph_pack::GraphPack &gp,
                               const std::string &hmm_set,
                               const std::filesystem::path &output_dir);
};

}

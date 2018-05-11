//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "domain_graph.hpp"

#include "common/pipeline/stage.hpp"

namespace hmmer {
class HMMMatcher;
class HMM;
struct hmmer_cfg;
}

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
private:
    std::vector<std::string> getFileVector(const std::string &hmm_files);
    void match_contigs_internal(hmmer::HMMMatcher &matcher, path_extend::BidirectionalPath* path,
                                const std::string &path_string, const hmmer::HMM &hmm, ContigAlnInfo &res, io::OFastaReadStream &oss_contig);
    ContigAlnInfo match_contigs(const path_extend::PathContainer &contig_paths,
                                const hmmer::HMM &hmm, const hmmer::hmmer_cfg &cfg,
                                const path_extend::ScaffoldSequenceMaker &scaffold_maker, io::OFastaReadStream &oss_contig);
};

}

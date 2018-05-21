//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "domain_matcher.hpp"

#include "hmm/hmmfile.hpp"
#include "hmm/hmmmatcher.hpp"

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"

#include "pipeline/graph_pack.hpp"

#include "sequence/aa.hpp"
#include "io/reads/osequencestream.hpp"

#include <string>
#include <vector>

namespace nrps {

static void match_contigs_internal(hmmer::HMMMatcher &matcher, path_extend::BidirectionalPath* path,
                                   const std::string &path_string,
                                   const std::string &type, const hmmer::HMM &hmm, ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    for (size_t shift = 0; shift < 3; ++shift) {
        std::string ref_shift = std::to_string(path->GetId()) + "_" + std::to_string(shift);
        std::string seq_aas = aa::translate(path_string.c_str() + shift);
        matcher.match(ref_shift.c_str(), seq_aas.c_str());
    }
    matcher.summarize();

    for (const auto &hit : matcher.hits()) {
        if (!hit.reported() || !hit.included())
            continue;

        for (const auto &domain : hit.domains()) {
            std::pair<int, int> seqpos = domain.seqpos();
            std::pair<int, int> seqpos2 = domain.hmmpos();
            INFO("First - " << seqpos2.first << ", second - " << seqpos2.second);
            INFO("First - " << seqpos.first << ", second - " << seqpos.second);
            int shift = hit.name()[strlen(hit.name()) - 1] - '0';
            seqpos.first = seqpos.first * 3  + shift;
            seqpos.second = seqpos.second * 3  + shift;
            std::string name(hit.name());
            oss_contig << io::SingleRead(name, path_string);
            INFO(name);
            INFO("First - " << seqpos.first << ", second - " << seqpos.second);
            res.push_back({type, name, unsigned(seqpos.first), unsigned(seqpos.second), path_string.substr(seqpos.first, seqpos.second - seqpos.first)});
        }
    }
}

static void match_contigs(const path_extend::PathContainer &contig_paths, const path_extend::ScaffoldSequenceMaker &scaffold_maker,
                          const std::string &type, const hmmer::HMM &hmm, const hmmer::hmmer_cfg &cfg,
                          ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    INFO("Total contigs: " << contig_paths.size());
    INFO("Model length - " << hmm.length());
    for (auto iter = contig_paths.begin(); iter != contig_paths.end(); ++iter) {
        hmmer::HMMMatcher matcher(hmm, cfg);
        path_extend::BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        std::string path_string = scaffold_maker.MakeSequence(*path);
        match_contigs_internal(matcher, path, path_string, type, hmm, res, oss_contig);

        path_extend::BidirectionalPath* conj_path = path->GetConjPath();
        if (conj_path->Length() <= 0)
            continue;
        std::string path_string_conj = scaffold_maker.MakeSequence(*conj_path);
        match_contigs_internal(matcher, conj_path, path_string_conj, type, hmm, res, oss_contig);
    }
}

ContigAlnInfo DomainMatcher::MatchDomains(debruijn_graph::conj_graph_pack &gp,
                                          const std::string &hmm_set,
                                          const std::string &output_dir) {
    if (fs::check_existence(output_dir + "/temp_anti/"))
        fs::remove_dir(output_dir + "/temp_anti/");
    fs::make_dirs(output_dir + "/temp_anti/");
    fs::make_dirs(output_dir + "/bgc_in_gfa/");

    ContigAlnInfo res;
    hmmer::hmmer_cfg hcfg;
    hcfg.E = 1.0e-9;
    hcfg.domE = 1.0e-9;
    std::vector<std::string> hmms;
    boost::split(hmms, hmm_set, boost::is_any_of(",;"), boost::token_compress_on);
    path_extend::ScaffoldSequenceMaker scaffold_maker(gp.g);
    path_extend::PathContainer broken_scaffolds;
    path_extend::ScaffoldBreaker(int(gp.g.k())).Break(gp.contig_paths, broken_scaffolds);

    io::OFastaReadStream oss_contig(output_dir + "/temp_anti/restricted_edges.fasta");
    for (const auto &file : hmms) {
        hmmer::HMMFile hmmfile(file);
        if (!hmmfile.valid())
        FATAL_ERROR("Error opening HMM file "<< file);
        auto hmmw = hmmfile.read();
        INFO("Matching contigs with " << file);

        std::string type = fs::filename(file);
        size_t dot = type.find_first_of(".");
        VERIFY(dot != std::string::npos);
        type = type.substr(0, dot);

        match_contigs(broken_scaffolds, scaffold_maker,
                      type, hmmw.get(), hcfg,
                      res, oss_contig);
    }

    return res;
}

}

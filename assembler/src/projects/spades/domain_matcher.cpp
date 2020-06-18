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
#include "utils/filesystem/path_helper.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>

namespace nrps {

static void match_contigs_internal(hmmer::HMMMatcher &matcher, path_extend::BidirectionalPath* path,
                                   const std::string &path_string,
                                   const std::string &type, ContigAlnInfo &res, io::OFastaReadStream &oss_contig, size_t model_length) {
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
            DEBUG("First - " << seqpos2.first << ", second - " << seqpos2.second);
            DEBUG("First - " << seqpos.first << ", second - " << seqpos.second);
            int shift = hit.name()[strlen(hit.name()) - 1] - '0';
            seqpos.first = seqpos.first * 3  + shift;
            seqpos.second = seqpos.second * 3  + shift;

            std::string name(hit.name());
#pragma omp critical
            {
                oss_contig << io::SingleRead(name, path_string);
            }
            DEBUG(name);
            DEBUG("First - " << seqpos.first << ", second - " << seqpos.second);
            res.push_back({type, name, unsigned(seqpos.first), unsigned(seqpos.second), path_string.substr(seqpos.first, seqpos.second - seqpos.first)});
        }
    }
    matcher.reset_top_hits();
}

static void match_contigs(const path_extend::PathContainer &contig_paths, const path_extend::ScaffoldSequenceMaker &scaffold_maker,
                          const std::string &type, const hmmer::HMM &hmm, const hmmer::hmmer_cfg &cfg,
                          ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    DEBUG("Total contigs: " << contig_paths.size());
    DEBUG("Model length - " << hmm.length());
    hmmer::HMMMatcher matcher(hmm, cfg);
    for (auto iter = contig_paths.begin(); iter != contig_paths.end(); ++iter) {
        path_extend::BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        std::string path_string = scaffold_maker.MakeSequence(*path);
        match_contigs_internal(matcher, path, path_string, type, res, oss_contig, hmm.length());

        path_extend::BidirectionalPath* conj_path = path->GetConjPath();
        if (conj_path->Length() <= 0)
            continue;
        std::string path_string_conj = scaffold_maker.MakeSequence(*conj_path);
        match_contigs_internal(matcher, conj_path, path_string_conj, type, res, oss_contig, hmm.length());
    }
}


static void ParseHMMFile(std::vector<hmmer::HMM> &hmms, const std::string &filename) {
    auto hmmfile = hmmer::open_file(filename);
    if (std::error_code ec = hmmfile.getError()) {
        FATAL_ERROR("Error opening HMM file "<< filename << ", reason: " << ec.message());
    }

    while (auto hmmw = hmmfile->read()) {
        if (std::error_code ec = hmmw.getError())
            FATAL_ERROR("Error reading HMM file "<< filename << ", reason: " << ec.message());
        hmms.emplace_back(std::move(hmmw.get()));
    }
}

ContigAlnInfo DomainMatcher::MatchDomains(debruijn_graph::GraphPack &gp,
                                          const std::string &hmm_set,
                                          const std::string &output_dir) {
    if (fs::check_existence(output_dir + "/temp_anti/"))
        fs::remove_dir(output_dir + "/temp_anti/");
    fs::make_dirs(output_dir + "/temp_anti/");
    fs::make_dirs(output_dir + "/bgc_in_gfa/");

    ContigAlnInfo res;
    hmmer::hmmer_cfg hcfg;
    hcfg.cut_ga = true;
    std::vector<std::string> hmm_files;
    boost::split(hmm_files, hmm_set, boost::is_any_of(",;"), boost::token_compress_on);
    path_extend::ScaffoldSequenceMaker scaffold_maker(gp.get<debruijn_graph::Graph>());
    path_extend::PathContainer broken_scaffolds;
    path_extend::ScaffoldBreaker(int(gp.k())).Break(gp.get<path_extend::PathContainer>(), broken_scaffolds);

    io::OFastaReadStream oss_contig(output_dir + "/temp_anti/restricted_edges.fasta");

    std::vector<hmmer::HMM> hmms;
    for (const auto &f : hmm_files)
        ParseHMMFile(hmms, f);
    
#   pragma omp parallel for
    for (size_t i = 0; i < hmms.size(); ++i) {
        ContigAlnInfo local_res;
        std::string name = hmms[i].name();
#       pragma omp critical
        {
            INFO("Matching contigs with " << name);
        }

        match_contigs(broken_scaffolds, scaffold_maker,
                      name, hmms[i], hcfg,
                      local_res, oss_contig);

#       pragma omp critical
        {
            INFO("Matches for '" << name << "': " << local_res.size());
            res.insert(res.end(), std::make_move_iterator(local_res.begin()), std::make_move_iterator(local_res.end()));
        }
    }

    INFO("Total domain matches: " << res.size());
    return res;
}

}

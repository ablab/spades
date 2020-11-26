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

extern "C" {
    #include "easel.h"
    #include "esl_sqio.h"
}

namespace nrps {

static void match_contigs_internal(hmmer::HMMMatcher &matcher, const path_extend::BidirectionalPath &path,
                                   const std::string &path_string,
                                   const std::string &type, const std::string &desc,
                                   ContigAlnInfo &res, io::OFastaReadStream &oss_contig, size_t model_length) {
    for (size_t shift = 0; shift < 3; ++shift) {
        std::string ref_shift = std::to_string(path.GetId()) + "_" + std::to_string(shift);
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
            if (seqpos2.second - seqpos2.first < model_length / 10) {
                DEBUG("Fragmented hit");
                continue;
            }
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
            res.push_back({name, type, desc,
                           unsigned(seqpos.first), unsigned(seqpos.second),
                           path_string.substr(seqpos.first, std::max(seqpos.second - seqpos.first, (int)path.g().k() + 1))});
        }
    }
    matcher.reset_top_hits();
}

static void match_contigs(const path_extend::PathContainer &contig_paths, const path_extend::ScaffoldSequenceMaker &scaffold_maker,
                          const hmmer::HMM &hmm, const hmmer::hmmer_cfg &cfg,
                          ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    DEBUG("Total contigs: " << contig_paths.size());
    DEBUG("Model length - " << hmm.length());
    hmmer::HMMMatcher matcher(hmm, cfg);
    for (auto iter = contig_paths.begin(); iter != contig_paths.end(); ++iter) {
        const path_extend::BidirectionalPath &path = iter.get();
        if (path.Length() <= 0)
            continue;

        std::string path_string = scaffold_maker.MakeSequence(path);
        match_contigs_internal(matcher, path, path_string,
                               hmm.name(), hmm.desc() ? hmm.desc() : "",
                               res, oss_contig, hmm.length());

        const path_extend::BidirectionalPath& conj_path = iter.getConjugate();
        if (conj_path.Length() <= 0)
            continue;

        std::string path_string_conj = scaffold_maker.MakeSequence(conj_path);
        match_contigs_internal(matcher, conj_path, path_string_conj,
                               hmm.name(), hmm.desc() ? hmm.desc() : "",
                               res, oss_contig, hmm.length());
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

static void ParseFASTAFile(std::vector<hmmer::HMM> &hmms, const std::string &filename) {
    hmmer::HMMSequenceBuilder builder(hmmer::Alphabet::AMINO, hmmer::ScoreSystem::Default);

    ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
    ESL_SQ         *qsq  = esl_sq_CreateDigital(abc);
    ESL_SQFILE     *qfp  = NULL;
    const char *qfile = filename.c_str();

    // Open the query sequence file in FASTA format
    int status = esl_sqfile_Open(qfile, eslSQFILE_FASTA, NULL, &qfp);
    if (status == eslENOTFOUND) {
        FATAL_ERROR("No such file " << filename);
    } else if (status == eslEFORMAT) {
        FATAL_ERROR("Format of " << filename << " unrecognized.");
    } else if (status == eslEINVAL) {
        FATAL_ERROR("Can't autodetect stdin or .gz.");
    } else if (status != eslOK) {
        FATAL_ERROR("Open of " << filename << " failed, code " << status);
    }

    // For each sequence, build a model and save it.
    while ((status = esl_sqio_Read(qfp, qsq)) == eslOK) {
        INFO("Converting " << qsq->name << ", len: " << qsq->n);
        hmms.push_back(builder.from_string(qsq));
        esl_sq_Reuse(qsq);
    }
    if (status != eslEOF) {
        FATAL_ERROR("Unexpected error " << status << " reading sequence file " << filename);
    }

    esl_sq_Destroy(qsq);
    esl_sqfile_Close(qfp);
    esl_alphabet_Destroy(abc);
}

ContigAlnInfo DomainMatcher::MatchDomains(debruijn_graph::GraphPack &gp,
                                          const std::string &hmm_set,
                                          const std::string &output_dir) {
    std::string tmp_dir = fs::append_path(output_dir, "temp_anti");
    if (fs::check_existence(tmp_dir))
        fs::remove_dir(tmp_dir);
    fs::make_dirs(tmp_dir);
    fs::make_dirs(fs::append_path(output_dir, "bgc_in_gfa"));

    ContigAlnInfo res;
    hmmer::hmmer_cfg hcfg;
    hcfg.cut_ga = true;
    std::vector<std::string> hmm_files;
    boost::split(hmm_files, hmm_set, boost::is_any_of(",;"), boost::token_compress_on);
    path_extend::ScaffoldSequenceMaker scaffold_maker(gp.get<debruijn_graph::Graph>());
    path_extend::PathContainer broken_scaffolds;
    path_extend::ScaffoldBreaker(int(gp.k())).Break(gp.get<path_extend::PathContainer>("exSPAnder paths"), broken_scaffolds);

    io::OFastaReadStream oss_contig(fs::append_path(output_dir, "restricted_edges.fasta"));

    std::vector<hmmer::HMM> hmms;
    for (const auto &f : hmm_files) {
        if (utils::ends_with(f, ".aa") || utils::ends_with(f, ".aa.gz")) {
            ParseFASTAFile(hmms, f);
            hcfg.E = hcfg.domE = 0.01;
        } else {
            ParseHMMFile(hmms, f);
        }
    }

    // Setup E-value search space size
    hcfg.Z = 3 * broken_scaffolds.size();

#   pragma omp parallel for
    for (size_t i = 0; i < hmms.size(); ++i) {
        ContigAlnInfo local_res;
        std::string name = hmms[i].name();
#       pragma omp critical
        {
            INFO("Matching contigs with " << name);
        }

        match_contigs(broken_scaffolds, scaffold_maker,
                      hmms[i], hcfg,
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

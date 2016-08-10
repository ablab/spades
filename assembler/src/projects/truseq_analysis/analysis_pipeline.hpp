//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"
#include <pipeline/stage.hpp>
#include "alignment_analyser.hpp"
#include "AlignmentAnalyserNew.hpp"

namespace spades {
    class VariationDetectionStage : public AssemblyStage {
    public:
        typedef debruijn_graph::config::debruijn_config::truseq_analysis Config;
    private:
        std::string output_file_;
        const Config &config_;
    public:
        VariationDetectionStage(string output_file, const Config &config);

        vector<io::SingleRead> ReadScaffolds(const string &scaffolds_file);

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *);

        DECL_LOGGER("AlignmntAnalysis")

        bool CheckEndVertex(debruijn_graph::DeBruijnGraph const &graph,
                                                             debruijn_graph::EdgeId id, size_t i);
    private:
        vector <alignment_analysis::ConsistentMapping> ExtractConsistentMappings(const vector<alignment_analysis::ConsistentMapping> &path);
    };

    void run_truseq_analysis();

}

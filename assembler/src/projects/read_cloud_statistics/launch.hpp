#pragma once

#include "common/pipeline/stage.hpp"
#include "projects/spades/repeat_resolving.hpp"
#include "read_cloud_statistics_extractor.hpp"
#include "scaffold_graph_construction_stage.hpp"
#include "scaffolder_analysis_stage.hpp"

namespace spades {
    void run_statistics_extractor() {
        INFO("Read cloud statistics extractor started");
        debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                                cfg::get().tmp_dir,
                                                cfg::get().ds.reads.lib_count(),
                                                cfg::get().ds.reference_genome,
                                                cfg::get().flanking_range,
                                                cfg::get().pos.max_mapping_gap,
                                                cfg::get().pos.max_gap_diff);
        conj_gp.kmer_mapper.Attach();
        StageManager statistics_extractor({cfg::get().developer_mode,
                              cfg::get().load_from,
                              cfg::get().output_saves});
        statistics_extractor.add(new debruijn_graph::RepeatResolution())
                .add(new debruijn_graph::ReadCloudStatisticsStage());
        INFO("Output directory: " << cfg::get().output_dir);
        statistics_extractor.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("Read cloud statistics extractor finished.");
    }

    void run_scaffolder_analysis() {
        INFO("Scaffolder analysis started");
        debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                                cfg::get().tmp_dir,
                                                cfg::get().ds.reads.lib_count(),
                                                cfg::get().ds.reference_genome,
                                                cfg::get().flanking_range,
                                                cfg::get().pos.max_mapping_gap,
                                                cfg::get().pos.max_gap_diff);
        conj_gp.kmer_mapper.Attach();
        StageManager statistics_extractor({cfg::get().developer_mode,
                                           cfg::get().load_from,
                                           cfg::get().output_saves});
        statistics_extractor.add(new debruijn_graph::RepeatResolution())
            .add(new debruijn_graph::ScaffoldGraphConstructionStage())
            .add(new debruijn_graph::ScaffolderAnalysisStage())
            .add(new debruijn_graph::ReadCloudStatisticsStage());
        INFO("Output directory: " << cfg::get().output_dir);
        statistics_extractor.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("Scaffolder analysis finished.");
    }
} //spades
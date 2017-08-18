#pragma once

#include "utils/logger/logger.hpp"
#include "common/pipeline/stage.hpp"
#include "old_extender_stats/general_barcode_statistics.hpp"

namespace debruijn_graph {
    class ReadCloudStatisticsStage : public spades::AssemblyStage {

    public:

        ReadCloudStatisticsStage() :
                AssemblyStage("Read cloud statistics extractor", "read_cloud_statistics_extractor") {
        }

        void save(const conj_graph_pack &, const std::string &, const char *) const { }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *);
        DECL_LOGGER("ReadCloudStatisticsStage")
    };

} //debruijn_graph

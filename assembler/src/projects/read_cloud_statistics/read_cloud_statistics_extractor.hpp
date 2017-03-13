#pragma once

#include "modules/path_extend/read_cloud_path_extend/extension_chooser_checker.hpp"
#include "utils/logger/logger.hpp"
#include "common/pipeline/stage.hpp"
#include "modules/path_extend/read_cloud_path_extend/general_barcode_statistics.hpp"

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

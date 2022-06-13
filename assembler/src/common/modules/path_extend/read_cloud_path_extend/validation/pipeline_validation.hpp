//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_construction_pipeline.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

class ScaffoldGraphPipelineValidator {
  public:
    typedef ScaffoldGraphConstructionPipeline Pipeline;
    ScaffoldGraphPipelineValidator(const string &path_to_reference,
                                   const ScaffoldingUniqueEdgeStorage &unique_storage,
                                   const graph_pack::GraphPack &gp);

    void ValidateStagesResults(const Pipeline &pipeline, const std::filesystem::path &output_path) const;
  private:
    const std::string path_to_reference_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
    const graph_pack::GraphPack &gp_;

    DECL_LOGGER("ScaffoldGraphPipelineValidator");
};

}
}
}
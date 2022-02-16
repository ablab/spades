//***************************************************************************
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "visualization.hpp"
#include "pipeline/genomic_quality.hpp"

namespace debruijn_graph {
template<class Graph>
class QualityEdgeLocalityPrintingRH : public QualityLoggingRemovalHandler<Graph> {
    typedef QualityLoggingRemovalHandler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    visualization::visualization_utils::LocalityPrintingRH<Graph> printing_rh_;
public:
    QualityEdgeLocalityPrintingRH(const Graph &g
            , const EdgeQuality<Graph> &quality_handler
            , const visualization::graph_labeler::GraphLabeler<Graph> &labeler
            , std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer
    , const std::string &output_folder, bool handle_all = false) :
    base(g, quality_handler, handle_all),
    printing_rh_(g, labeler, colorer, output_folder)
    {}

    void HandlePositiveQuality(EdgeId e) override {
        printing_rh_.HandleDelete(e, "_" + std::to_string(this->quality_handler().quality(e)));
    }

private:
    DECL_LOGGER("QualityEdgeLocalityPrintingRH");
};

}

//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//TODO: move into 'pipeline'?

#pragma once

#include <memory>

#include "graph.hpp"
#include "coverage.hpp"
#include "edge_index.hpp"
#include "kmer_mapper.hpp"
#include "paired_index.hpp"
#include "positions.hpp"
#include "pipeline/graph_pack.hpp"

namespace io {

//TODO: get rid of ad-hoc component processing

template<typename Graph>
class GraphPackIO : public IOBase<debruijn_graph::graph_pack<Graph>> {
public:
    typedef typename debruijn_graph::graph_pack<Graph> Type;
    void Save(const std::string &basename, const Type &gp) override {
        //Save basic graph
        graph_io_.Save(basename, gp.g);
        const auto &mapper = graph_io_.GetEdgeMapper(); //TODO: get rid of this unused parameter

        //Save coverage
        CoverageIO<CoverageIndex<Graph>>(mapper)
                .Save(basename, gp.g.coverage_index());

        if (gp.edge_pos.IsAttached()) { //Save edge positions
            EdgePositionsIO<Graph>(mapper)
                    .Save(basename, gp.edge_pos);
        }

        if (gp.index.IsAttached()) { //Save kmer edge index
            EdgeIndexIO<Graph>()
                    .Save(basename, gp.index);
        }

        if (gp.kmer_mapper.IsAttached()) { //Save kmer mapper
            KmerMapperIO<Graph>()
                    .Save(basename, gp.kmer_mapper);
        }

        if (gp.flanking_cov.IsAttached()) { //Save flanking coverage
            CoverageIO<FlankingCoverage<Graph>>(mapper)
                    .Save(basename, gp.flanking_cov);
        }
    }

    void Load(const std::string &basename, Type &gp) override {
        //Load basic graph
        graph_io_.Load(basename, gp.g);
        const auto &mapper = graph_io_.GetEdgeMapper();

        //Load coverage
        CoverageIO<CoverageIndex<Graph>>(mapper)
                .Load(basename, gp.g.coverage_index());

        //Load edge positions
        EdgePositionsIO<Graph>(mapper)
                .Load(basename, gp.edge_pos);

        //Load kmer edge index
        EdgeIndexIO<Graph>()
                .Load(basename, gp.index);

        if (gp.kmer_mapper.IsAttached()) { //Load kmer mapper
            KmerMapperIO<Graph>()
                    .Load(basename, gp.kmer_mapper);
        }

        //Load flanking coverage
        CoverageIO<FlankingCoverage<Graph>>(mapper)
                .Load(basename, gp.flanking_cov);
    }

protected:
    GraphIO<Graph> graph_io_;
};

template<typename Graph>
class FullPackIO : public GraphPackIO<Graph> {
public:
    typedef GraphPackIO<Graph> base;
    typedef typename debruijn_graph::graph_pack<Graph> Type;
    void Save(const std::string &basename, const Type &gp) override {
        //Save basic graph
        base::Save(basename, gp);
        const auto &mapper = this->graph_io_.GetEdgeMapper(); //TODO: get rid of this unused parameter

        //Save unclustered paired indices
        using namespace omnigraph::de;
        PairedIndicesIO<UnclusteredPairedInfoIndexT<Graph>>(mapper)
                .Save(basename, gp.paired_indices);

        { //Save clustered & scaffolding indices
            PairedIndicesIO<PairedInfoIndexT<Graph>> io(mapper);
            io.Save(basename + "_cl", gp.clustered_indices);
            io.Save(basename + "_scf", gp.scaffolding_indices);
        }

        { //Save long reads
        }

        gp.ginfo.Save(basename + ".ginfo");
    }

    void Load(const std::string &basename, Type &gp) override {
        //Load basic graph
        base::Load(basename, gp);
        const auto &mapper = this->graph_io_.GetEdgeMapper();

        //Load paired indices
        using namespace omnigraph::de;
        PairedIndicesIO<UnclusteredPairedInfoIndexT<Graph>>(mapper)
                .Load(basename, gp.paired_indices);

        { //Load clustered & scaffolding indices
            PairedIndicesIO<PairedInfoIndexT<Graph>> io(mapper);
            io.Load(basename + "_cl", gp.clustered_indices);
            io.Load(basename + "_scf", gp.scaffolding_indices);
        }

        { //Load long reads
        }

        gp.ginfo.Load(basename + ".ginfo");
    }
};

}

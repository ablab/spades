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
#include "long_reads.hpp"
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

        //Save coverage
        CoverageIO<Graph>()
                .Save(basename, gp.g.coverage_index());

        if (gp.edge_pos.IsAttached()) { //Save edge positions
            EdgePositionsIO<Graph>()
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
            FlankingCoverageIO<Graph>()
                    .Save(basename, gp.flanking_cov);
        }
    }

    bool Load(const std::string &basename, Type &gp) override {
        //Load basic graph
        bool loaded = graph_io_.Load(basename, gp.g);
        VERIFY(loaded);
        const auto &mapper = graph_io_.GetEdgeMapper();

        //Load coverage
        loaded = CoverageIO<Graph>()
                .Load(basename, gp.g.coverage_index(), mapper);
        VERIFY(loaded);

        //Load edge positions
        VERIFY(!gp.edge_pos.IsAttached());
        gp.edge_pos.Attach();
        if (!EdgePositionsIO<Graph>()
                .Load(basename, gp.edge_pos, mapper)) {
            INFO("No saved positions");
        }

        //Load kmer edge index
        if (!EdgeIndexIO<Graph>()
                .Load(basename, gp.index)) {
            WARN("Cannot load edge index, kmer coverages will be missed");
            gp.index.Refill();
        }

        if (gp.kmer_mapper.IsAttached()) { //Load kmer mapper
            if (!KmerMapperIO<Graph>()
                    .Load(basename, gp.kmer_mapper)) {
                WARN("Cannot load kmer_mapper, information on projected kmers will be missed");
            }
        }

        //Load flanking coverage
        if (!FlankingCoverageIO<Graph>()
                .Load(basename, gp.flanking_cov, mapper)) {
            WARN("Cannot load flanking coverage, flanking coverage will be recovered from index");
            gp.flanking_cov.Fill(gp.index.inner_index());
        }

        return true;
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

        //Save unclustered paired indices
        using namespace omnigraph::de;
        PairedIndicesIO<UnclusteredPairedInfoIndexT<Graph>>()
                .Save(basename, gp.paired_indices);

        { //Save clustered & scaffolding indices
            PairedIndicesIO<PairedInfoIndexT<Graph>> io;
            io.Save(basename + "_cl", gp.clustered_indices);
            io.Save(basename + "_scf", gp.scaffolding_indices);
        }

        //Save long reads
        LongReadsIO<Graph>().Save(basename, gp.single_long_reads);

        gp.ginfo.Save(basename + ".ginfo");
    }

    bool Load(const std::string &basename, Type &gp) override {
        //Load basic graph
        bool loaded = base::Load(basename, gp);
        VERIFY(loaded);
        const auto &mapper = this->graph_io_.GetEdgeMapper();

        //Load paired indices
        using namespace omnigraph::de;
        PairedIndicesIO<UnclusteredPairedInfoIndexT<Graph>>()
                .Load(basename, gp.paired_indices, mapper);

        { //Load clustered & scaffolding indices
            PairedIndicesIO<PairedInfoIndexT<Graph>> io;
            io.Load(basename + "_cl", gp.clustered_indices, mapper);
            io.Load(basename + "_scf", gp.scaffolding_indices, mapper);
        }

        //Load long reads
        LongReadsIO<Graph>().Load(basename, gp.single_long_reads, mapper);

        gp.ginfo.Load(basename + ".ginfo");

        return true;
    }
};

}

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

namespace binary {

//TODO: get rid of ad-hoc component processing

template<typename Graph>
class GraphPackIO : public IOBase<debruijn_graph::graph_pack<Graph>> {
public:
    typedef typename debruijn_graph::graph_pack<Graph> Type;

    void Save(const std::string &basename, const Type &gp) override {
        //1. Save basic graph
        graph_io_.Save(basename, gp.g);

        //2. Save coverage
        CoverageIO<Graph>().Save(basename, gp.g.coverage_index());

        //3. Save edge positions
        SaveAttached(basename, gp.edge_pos);

        //4. Save kmer edge index
        SaveAttached(basename, gp.index);

        //5. Save kmer mapper
        SaveAttached(basename, gp.kmer_mapper);

        //6. Save flanking coverage
        SaveAttached(basename, gp.flanking_cov);
    }

    bool Load(const std::string &basename, Type &gp) override {
        //1. Load basic graph
        bool loaded = graph_io_.Load(basename, gp.g);
        VERIFY(loaded);
        const auto &mapper = graph_io_.GetEdgeMapper();

        //2. Load coverage
        loaded = CoverageIO<Graph>().Load(basename, gp.g.coverage_index(), mapper);
        VERIFY(loaded);

        //3. Load edge positions
        LoadAttached(basename, gp.edge_pos, mapper);

        //4. Load kmer edge index
        LoadAttached(basename, gp.index);

        //5. Load kmer mapper
        LoadAttached(basename, gp.kmer_mapper);

        //6. Load flanking coverage
        LoadAttached(basename, gp.flanking_cov, mapper);

        return true;
    }

protected:
    GraphIO<Graph> graph_io_;

    template<typename T>
    void SaveAttached(const std::string &basename, const T &component) {
        if (component.IsAttached()) {
            typename IOTraits<T>::Type io;
            io.Save(basename, component);
        }
    }

    template<typename T, typename... Env>
    void LoadAttached(const std::string &basename, T &component, const Env &... env) {
        if (component.IsAttached())
            component.Detach();
        typename IOTraits<T>::Type io;
        if (io.Load(basename, component, env...))
            component.Attach();
    }
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

} // namespace binary

} // namespace io

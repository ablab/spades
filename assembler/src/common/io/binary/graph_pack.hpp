//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//TODO: move into 'pipeline'?

#pragma once

#include <memory>

#include "basic.hpp"
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
class BasePackIO : public IOBase<debruijn_graph::graph_pack<Graph>> {
public:
    typedef typename debruijn_graph::graph_pack<Graph> Type;

    void Save(const std::string &basename, const Type &gp) override {
        //1. Save basic graph
        graph_io_.Save(basename, gp.g);

        //2. Save edge positions
        SaveAttached(basename, gp.edge_pos);

        //3. Save kmer edge index
        SaveAttached(basename, gp.index);

        //4. Save kmer mapper
        SaveAttached(basename, gp.kmer_mapper);

        //5. Save flanking coverage
        SaveAttached(basename, gp.flanking_cov);
    }

    bool Load(const std::string &basename, Type &gp) override {
        //1. Load basic graph
        graph_io_.Load(basename, gp.g);

        //2. Load edge positions
        LoadAttached(basename, gp.edge_pos);

        //3. Load kmer edge index
        LoadAttached(basename, gp.index);

        //4. Load kmer mapper
        LoadAttached(basename, gp.kmer_mapper);

        //5. Load flanking coverage
        LoadAttached(basename, gp.flanking_cov);

        return true;
    }

protected:
    BasicGraphIO<Graph> graph_io_;

    template<typename T>
    void SaveAttached(const std::string &basename, const T &component) {
        if (component.IsAttached()) {
            typename IOTraits<T>::Type io;
            io.Save(basename, component);
        }
    }

    template<typename T>
    void LoadAttached(const std::string &basename, T &component) {
        if (component.IsAttached())
            component.Detach();
        typename IOTraits<T>::Type io;
        if (io.Load(basename, component))
            component.Attach();
    }
};

template<typename Graph>
class FullPackIO : public BasePackIO<Graph> {
public:
    typedef BasePackIO<Graph> base;
    typedef typename debruijn_graph::graph_pack<Graph> Type;
    void Save(const std::string &basename, const Type &gp) override {
        //1. Save basic graph
        base::Save(basename, gp);

        //2. Save unclustered paired indices
        using namespace omnigraph::de;
        io::binary::Save(basename, gp.paired_indices);

        //3. Save clustered & scaffolding indices
        io::binary::Save(basename + "_cl", gp.clustered_indices);
        io::binary::Save(basename + "_scf", gp.scaffolding_indices);

        //4. Save long reads
        io::binary::Save(basename, gp.single_long_reads);

        //5. Save genome info
        gp.ginfo.Save(basename + ".ginfo");
    }

    bool Load(const std::string &basename, Type &gp) override {
        //1. Load basic graph
        bool loaded = base::Load(basename, gp);
        VERIFY(loaded);

        //2. Load paired indices
        using namespace omnigraph::de;
        io::binary::Load(basename, gp.paired_indices);

        //3. Load clustered & scaffolding indices
        io::binary::Load(basename + "_cl", gp.clustered_indices);
        io::binary::Load(basename + "_scf", gp.scaffolding_indices);

        //4. Load long reads
        io::binary::Load(basename, gp.single_long_reads);

        //5. Load genome info
        gp.ginfo.Load(basename + ".ginfo");

        return true;
    }
};

} // namespace binary

} // namespace io

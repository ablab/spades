//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_pack.hpp"

#include "basic.hpp"
#include "coverage.hpp"
#include "edge_index.hpp"
#include "genomic_info.hpp"
#include "kmer_mapper.hpp"
#include "long_reads.hpp"
#include "ss_coverage.hpp"
#include "paired_index.hpp"
#include "positions.hpp"
#include "trusted_paths.hpp"

namespace io {

namespace binary {

namespace {

class Saver {
    const std::string &basename;
    const BasePackIO::Type &gp;
    std::ofstream infoStream;
public:
    Saver(const std::string &basename, const BasePackIO::Type &gp)
        : basename(basename)
        , gp(gp)
        , infoStream(basename + ".att")
    {}

    /**
     * @brief  Saves the component only if it was attached.
     *         Also adds its attachment flag to the attached metadata.
     */
    template<class T>
    void Save() {
        const auto &component = gp.get<T>();
        io::binary::BinWrite<char>(infoStream, component.IsAttached());
        if (component.IsAttached()) {
            typename IOTraits<T>::Type io;
            io.Save(basename, component); 
        }
    }
};

class BinWriter {
    std::ostream &os;
    const BasePackIO::Type &gp;
public:
    BinWriter(std::ostream &os, const BasePackIO::Type &gp)
        : os(os)
        , gp(gp)
    {}

    /**
     * @brief  Write the component only if it was attached.
     *         Also adds its attachment flag to the attached metadata.
     */
    template<class T>
    void Write() {
        const auto &component = gp.get<T>();
        io::binary::BinWrite<char>(os, component.IsAttached());
        if (component.IsAttached()) {
            typename IOTraits<T>::Type io;
            io.BinWrite(os, component); 
        }
    }
};

class Loader {
    const std::string &basename;
    BasePackIO::Type &gp;
    std::ifstream infoStream;
public:
    Loader(const std::string &basename, BasePackIO::Type &gp)
        : basename(basename)
        , gp(gp)
        , infoStream(fs::open_file(basename + ".att", std::ios::binary))
    {}

    /**
     * @brief  Restores the attachment flag of the component. Then loads it only if it was attached.
     */
    template<class T>
    void Load() {
        INFO("Trying to load " << typeid(T).name());
        if (!io::binary::BinRead<char>(infoStream)) {
            INFO("Not attached, skipping");
            return;
        }
        auto &component = gp.get_mutable<T>();
        if (component.IsAttached())
            component.Detach();
        typename IOTraits<T>::Type io;
        bool loaded = io.Load(basename, component);
        VERIFY(loaded);
        component.Attach();
    }
};

class BinReader {
    std::istream &is;
    BasePackIO::Type &gp;
public:
    BinReader(std::istream &is, BasePackIO::Type &gp)
        : is(is)
        , gp(gp)
    {}

    /**
     * @brief  Restores the attachment flag of the component. Then reads it only if it was attached.
     */
    template<class T>
    void Read() {
        INFO("Trying to read " << typeid(T).name());
        if (!io::binary::BinRead<char>(is)) {
            INFO("Not attached, skipping");
            return;
        }
        auto &component = gp.get_mutable<T>();
        if (component.IsAttached())
            component.Detach();
        typename IOTraits<T>::Type io;
        auto loaded = io.BinRead(is, component);
        VERIFY(loaded);
        component.Attach();
    }
};

/**
 * @brief  Saves the component.
 */
template<typename T>
void SaveComponent(const std::string &basename, const BasePackIO::Type &gp, const std::string &name = "") {
    const auto &component = gp.get<T>(name);
    io::binary::Save(basename, component);
}

/**
 * @brief  Loads an arbitrary component.
 */
template<typename T>
void LoadComponent(const std::string &basename, BasePackIO::Type &gp, const std::string &name = "") {
    auto &component = gp.get_mutable<T>(name);
    io::binary::Load(basename, component);
}

/**
 * @brief  Writes the component in binary mode.
 */
template<typename T>
void BinWriteComponent(std::ostream &os, const BasePackIO::Type &gp, const std::string &name = "") {
    const auto &component = gp.get<T>(name);
    io::binary::Write(os, component);
}

/**
 * @brief  Reads an arbitrary component in binary mode.
 */
template<typename T>
void BinReadComponent(std::istream &is, BasePackIO::Type &gp, const std::string &name = "") {
    auto &component = gp.get_mutable<T>(name);
    io::binary::Read(is, component);
}

} // namespace

void BasePackIO::Save(const std::string &basename, const Type &gp) {
    Saver saver(basename, gp);

    using namespace omnigraph;
    using namespace debruijn_graph;

    if (gp.invalidated<Graph>()) {
        //1. Save basic graph with coverage
        const auto &g = gp.get<Graph>();
        graph_io_.Save(basename, g);
    }

    //2. Save edge positions
    saver.Save<EdgesPositionHandler<Graph>>();

    //3. Save kmer edge index
    saver.Save<EdgeIndex<Graph>>();

    //4. Save kmer mapper
    saver.Save<KmerMapper<Graph>>();

    //5. Save flanking coverage
    saver.Save<FlankingCoverage<Graph>>();
}

bool BasePackIO::Load(const std::string &basename, Type &gp) {
    Loader loader(basename, gp);

    using namespace omnigraph;
    using namespace debruijn_graph;

    //1. Load basic graph with coverage
    auto &g = gp.get_mutable<Graph>();
    graph_io_.Load(basename, g);

    //2. Load edge positions
    loader.Load<EdgesPositionHandler<Graph>>();

    //3. Load kmer edge index
    loader.Load<EdgeIndex<Graph>>();

    //4. Load kmer mapper
    loader.Load<KmerMapper<Graph>>();

    //5. Load flanking coverage
    loader.Load<FlankingCoverage<Graph>>();

    return true;
}

void BasePackIO::BinWrite(std::ostream &os, const Type &gp)  {
    BinWriter writer(os, gp);

    using namespace omnigraph;
    using namespace debruijn_graph;

    //1. Write basic graph
    graph_io_.BinWrite(os, gp.get<Graph>());

    //2. Write edge positions
    writer.Write<EdgesPositionHandler<Graph>>();

    //3. Write kmer edge index
    writer.Write<EdgeIndex<Graph>>();

    //4. Write kmer mapper
    writer.Write<KmerMapper<Graph>>();

    //5. Write flanking coverage
    writer.Write<FlankingCoverage<Graph>>();
}

bool BasePackIO::BinRead(std::istream &is, Type &gp) {
    BinReader reader(is, gp);

    using namespace omnigraph;
    using namespace debruijn_graph;

    //1. Read basic graph
    graph_io_.BinRead(is, gp.get_mutable<Graph>());

    //2. Read edge positions
    reader.Read<EdgesPositionHandler<Graph>>();

    //3. Read kmer edge index
    reader.Read<EdgeIndex<Graph>>();

    //4. Read kmer mapper
    reader.Read<KmerMapper<Graph>>();

    //5. Read flanking coverage
    reader.Read<FlankingCoverage<Graph>>();

    return true;
}

void FullPackIO::Save(const std::string &basename, const Type &gp) {
    using namespace omnigraph::de;
    using namespace debruijn_graph;

    //1. Save basic graph pack
    base::Save(basename, gp);

    //2. Save unclustered paired indices
    SaveComponent<UnclusteredPairedInfoIndicesT<Graph>>(basename, gp);

    //3. Save clustered indices
    SaveComponent<PairedInfoIndicesT<Graph>>(basename + "_cl", gp, "clustered_indices");

    //4. Save scaffolding indices
    SaveComponent<PairedInfoIndicesT<Graph>>(basename + "_scf", gp, "scaffolding_indices");

    //5. Save long reads
    SaveComponent<LongReadContainer<Graph>>(basename, gp);

    //6. Save genomic info
    SaveComponent<GenomicInfo>(basename, gp);

    //7. Save SS coverage
    SaveComponent<SSCoverageContainer>(basename, gp);

    //8. Save trusted paths
    SaveComponent<path_extend::TrustedPathsContainer>(basename, gp);
}

bool FullPackIO::Load(const std::string &basename, Type &gp) {
    using namespace omnigraph::de;
    using namespace debruijn_graph;

    //1. Load basic graph pack
    base::Load(basename, gp);

    //2. Load paired indices
    using namespace omnigraph::de;
    LoadComponent<UnclusteredPairedInfoIndicesT<Graph>>(basename, gp);

    //3. Load clustered indices
    LoadComponent<PairedInfoIndicesT<Graph>>(basename + "_cl", gp, "clustered_indices");

    //4. Load scaffolding indices
    LoadComponent<PairedInfoIndicesT<Graph>>(basename + "_scf", gp, "scaffolding_indices");

    //5. Load long reads
    LoadComponent<LongReadContainer<Graph>>(basename, gp);

    //6. Load genomic info
    LoadComponent<GenomicInfo>(basename, gp);

    //7. Load SS coverage
    LoadComponent<SSCoverageContainer>(basename, gp);

    //8. Load trusted paths
    LoadComponent<path_extend::TrustedPathsContainer>(basename, gp);

    return true;
}

void FullPackIO::BinWrite(std::ostream &os, const Type &gp) {
    using namespace omnigraph::de;
    using namespace debruijn_graph;

    //1. Write basic graph
    base::BinWrite(os, gp);

    //2. Write unclustered paired indices
    using namespace omnigraph::de;
    BinWriteComponent<UnclusteredPairedInfoIndicesT<Graph>>(os, gp);

    //3. Write clustered indices
    BinWriteComponent<PairedInfoIndicesT<Graph>>(os, gp, "clustered_indices");

    //4. Write scaffolding indices
    BinWriteComponent<PairedInfoIndicesT<Graph>>(os, gp, "scaffolding_indices");

    //5. Write long reads
    BinWriteComponent<LongReadContainer<Graph>>(os, gp);

    //6. Write genomic info
    BinWriteComponent<GenomicInfo>(os, gp);

    //7. Write SS coverage
    BinWriteComponent<SSCoverageContainer>(os, gp);
}

bool FullPackIO::BinRead(std::istream &is, Type &gp) {
    using namespace omnigraph::de;
    using namespace debruijn_graph;

    //1. Read basic graph
    bool loaded = base::BinRead(is, gp);
    VERIFY(loaded);

    //2. Read unclustered paired indices
    using namespace omnigraph::de;
    BinReadComponent<UnclusteredPairedInfoIndicesT<Graph>>(is, gp);

    //3. Read clustered indices
    BinReadComponent<PairedInfoIndicesT<Graph>>(is, gp, "clustered_indices");

    //4. Read scaffolding indices
    BinReadComponent<PairedInfoIndicesT<Graph>>(is, gp, "scaffolding_indices");

    //5. Read long reads
    BinReadComponent<LongReadContainer<Graph>>(is, gp);

    //6. Read genomic info
    BinReadComponent<GenomicInfo>(is, gp);

    //7. Read SS coverage
    BinReadComponent<SSCoverageContainer>(is, gp);

    return true;
}

} // namespace binary

} // namespace io

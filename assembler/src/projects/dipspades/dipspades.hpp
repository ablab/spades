//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "io/reads/splitting_wrapper.hpp"
#include "modules/graph_construction.hpp"
#include "pipeline/stage.hpp"

#include "dipspades_config.hpp"

#include "polymorphic_bulge_remover/polymorphic_bulge_remover.hpp"
#include "consensus_contigs_constructor/consensus_contigs_constructor.hpp"
#include "haplotype_assembly/haplotype_assembler.hpp"
#include "kmer_gluing/equal_sequence_gluer.hpp"

using namespace debruijn_graph;
using namespace spades;

namespace dipspades {

void construct_graph_from_contigs(debruijn_graph::conj_graph_pack &graph_pack){
    auto fnames = GetAllLinesFromFile(dsp_cfg::get().io.haplocontigs);
    ReadStreamList<SingleRead> streams;
    for(auto fname = fnames.begin(); fname != fnames.end(); fname++)
        if(fname_valid(*fname)){
            INFO("Addition of contigs from " << *fname);
            streams.push_back(io::SplittingWrap(EasyStream(*fname, true)));
        }

    INFO("Construction of the de Bruijn graph with K=" << dsp_cfg::get().bp.K);
    config::debruijn_config::construction params;
    params.con_mode = config::construction_mode::extention;
    params.early_tc.enable = false;
    params.early_tc.length_bound = 10;
    params.keep_perfect_loops = true;
    params.read_buffer_size = dsp_cfg::get().bp.read_buffer_size;

    ConstructGraphWithCoverage(params,
            streams,
            graph_pack.g,
            graph_pack.index,
            graph_pack.flanking_cov);
}


class DipSPAdesStorage{
public:
    BaseHistogram<size_t> bulge_len_histogram;
    ContigStoragePtr default_storage;
    ContigStoragePtr composite_storage;
    CorrectionResult redundancy_map;
};


class DipSPAdes : public CompositeStage<DipSPAdesStorage> {
    DipSPAdesStorage dsp_params_;
public:
    DipSPAdes() : CompositeStage<DipSPAdesStorage>("dipSPAdes", "dipspades") { }

    void load(debruijn_graph::conj_graph_pack&,
            const std::string &,
            const char*) { }

    void save(const debruijn_graph::conj_graph_pack&,
            const std::string &,
            const char*) const { }

    virtual ~DipSPAdes() { }
};

class ContigGraphConstructionStage : public DipSPAdes::Phase {
public:
    ContigGraphConstructionStage() :
        DipSPAdes::Phase("Construction of graph from contigs", "contig_graph_construction") { }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char*) {
        construct_graph_from_contigs(graph_pack);
    }

    void load(debruijn_graph::conj_graph_pack& gp,
            const std::string &load_from,
            const char* prefix) {
        std::string p = path::append_path(load_from, prefix == NULL ? id() : prefix);
        INFO("Loading current state from " << p);
        debruijn_graph::graphio::ScanAll(p, gp, false);

    }

    void save(const debruijn_graph::conj_graph_pack& gp,
            const std::string & save_to,
            const char* prefix) const {
        std::string p = path::append_path(save_to, prefix == NULL ? id() : prefix);
        INFO("Saving current state to " << p);
        debruijn_graph::graphio::PrintAll(p, gp);
    }

    virtual ~ContigGraphConstructionStage() { }
};

class PolymorphicBulgeRemoverStage : public DipSPAdes::Phase {
public:
    PolymorphicBulgeRemoverStage() :
        DipSPAdes::Phase("Polymorphic bulge remover", "polymorphic_br") { }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char*){
        if(dsp_cfg::get().pbr.enabled){
            PolymorphicBulgeRemover(graph_pack, storage().bulge_len_histogram).Run();
              INFO("Consensus graph was constructed");
        }
    }

    void load(debruijn_graph::conj_graph_pack& gp,
            const std::string &load_from,
            const char* prefix) {
        std::string p = path::append_path(load_from, prefix == NULL ? id() : prefix);
        INFO("Loading current state from " << p);
        debruijn_graph::graphio::ScanAll(p, gp, false);
        INFO("Loading histogram of bulge length");
        INFO("loading from " << p + ".hist");
        storage().bulge_len_histogram.LoadFrom(p + ".hist");
    }

    void save(const debruijn_graph::conj_graph_pack& gp,
            const std::string & save_to,
            const char* prefix) const {
        std::string p = path::append_path(save_to, prefix == NULL ? id() : prefix);
        INFO("Saving current state to " << p);
        debruijn_graph::graphio::PrintAll(p, gp);
        storage().bulge_len_histogram.SaveToFile(p + ".hist");
    }

    virtual ~PolymorphicBulgeRemoverStage() { }
};

class EqualKmerGluingStage : public DipSPAdes::Phase {
public:
    EqualKmerGluingStage() :
            DipSPAdes::Phase("Equal k-mer gluing", "kmer_gluer") { }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char*) {
        INFO("Glueing equal kmers starts");
        EqualSequencesGluer<Graph>(graph_pack.g, graph_pack.index).GlueEqualKmers();
        INFO("Glueing equal kmers ends");
    }

    void load(debruijn_graph::conj_graph_pack& gp,
              const std::string &load_from,
              const char* prefix) {
        std::string p = path::append_path(load_from, prefix == NULL ? id() : prefix);
        INFO("Loading current state from " << p);
        debruijn_graph::graphio::ScanAll(p, gp, false);
        INFO("Loading histogram of bulge length");
        INFO("loading from " << p + ".hist");
        storage().bulge_len_histogram.LoadFrom(p + ".hist");
    }

    void save(const debruijn_graph::conj_graph_pack& gp,
              const std::string & save_to,
              const char* prefix) const {
        std::string p = path::append_path(save_to, prefix == NULL ? id() : prefix);
        INFO("Saving current state to " << p);
        debruijn_graph::graphio::PrintAll(p, gp);
        storage().bulge_len_histogram.SaveToFile(p + ".hist");
    }

    virtual ~EqualKmerGluingStage() { }
};

class ConsensusConstructionStage : public DipSPAdes::Phase {
public:
    ConsensusConstructionStage() :
        DipSPAdes::Phase("Consensus contigs construction", "consensus_construction") { }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char*){
        if(dsp_cfg::get().cc.enabled){
            ConsensusContigsConstructor consensus_constructor(graph_pack, storage().bulge_len_histogram);
            consensus_constructor.Run();
            storage().composite_storage = consensus_constructor.CompositeContigsStorage();
            storage().default_storage = consensus_constructor.DefaultContigsStorage();
            storage().redundancy_map = consensus_constructor.RedundancyResult();
        }
    }

    void load(debruijn_graph::conj_graph_pack& gp,
            const std::string &load_from,
            const char* prefix) {
        std::string p = path::append_path(load_from, prefix == NULL ? id() : prefix);
        INFO("Loading current state from " << p);
        debruijn_graph::graphio::ScanAll(p, gp, false);
    }

    void save(const debruijn_graph::conj_graph_pack& gp,
            const std::string & save_to,
            const char* prefix) const {
        std::string p = path::append_path(save_to, prefix == NULL ? id() : prefix);
        INFO("Saving current state to " << p);
        debruijn_graph::graphio::PrintAll(p, gp);
        storage().bulge_len_histogram.SaveToFile(p + ".hist");
    }

    virtual ~ConsensusConstructionStage() { }
};

class HaplotypeAssemblyStage : public DipSPAdes::Phase {
public:
    HaplotypeAssemblyStage() :
        DipSPAdes::Phase("Haplotype assembly", "haplotype_assembly") { }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char*) {
        if(!storage().composite_storage || !storage().default_storage)
            return;
        if(storage().composite_storage->Size() == 0 || storage().default_storage->Size() == 0)
            return;
        INFO("Diploid graph construction");
        conj_graph_pack double_graph_pack(graph_pack.k_value, dsp_cfg::get().io.tmp_dir,
                                          dsp_cfg::get().io.num_libraries, "");
        construct_graph_from_contigs(double_graph_pack);
        HaplotypeAssembler(graph_pack, double_graph_pack, storage().default_storage,
                           storage().composite_storage, storage().redundancy_map).Run();
    }

    void load(debruijn_graph::conj_graph_pack&,
            const std::string &,
            const char*) { }

    void save(const debruijn_graph::conj_graph_pack&,
            const std::string &,
            const char*) const { }

    virtual ~HaplotypeAssemblyStage() { }
};
void run_dipspades() {
    INFO("dipSPAdes started");

    debruijn_graph::conj_graph_pack conj_gp(
            dsp_cfg::get().bp.K,
            dsp_cfg::get().io.tmp_dir,
            dsp_cfg::get().io.num_libraries,
            "", // reference genome
            1); // flanking range

    conj_gp.kmer_mapper.Attach();

    StageManager DS_Manager ( {dsp_cfg::get().rp.developer_mode,
                               dsp_cfg::get().io.saves,
                               dsp_cfg::get().io.output_saves} );
    auto ds_phase = new DipSPAdes();
    ds_phase -> add(new ContigGraphConstructionStage()) ->
            add(new PolymorphicBulgeRemoverStage()) ->
            add(new EqualKmerGluingStage()) ->
            add(new ConsensusConstructionStage());
    if(dsp_cfg::get().ha.ha_enabled) {
        ds_phase->add(new HaplotypeAssemblyStage());
    }

    DS_Manager.add(ds_phase);
    DS_Manager.run(conj_gp, dsp_cfg::get().rp.entry_point.c_str());
    INFO("dipSPAdes finished");
}

}

#pragma once

#include "pipeline/stage.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "algorithms/simplification/tip_clipper.hpp"
#include "projects/mts/contig_abundance.hpp"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

namespace debruijn_graph {

struct SeriesAnalysisConfig {
    uint k, sample_cnt;
    std::string kmer_mult, bin, bin_prof;
};

}

namespace llvm { namespace yaml {

template<> struct MappingTraits<debruijn_graph::SeriesAnalysisConfig> {
    static void mapping(IO& io, debruijn_graph::SeriesAnalysisConfig& cfg) {
        io.mapRequired("k", cfg.k);
        io.mapRequired("sample_cnt", cfg.sample_cnt);
        io.mapRequired("kmer_mult", cfg.kmer_mult);
        io.mapRequired("bin", cfg.bin);
        io.mapRequired("bin_prof", cfg.bin_prof);
    }
};

} }

namespace debruijn_graph {

template<class graph_pack>
shared_ptr<omnigraph::visualization::GraphColorer<typename graph_pack::graph_t>> DefaultGPColorer(
    const graph_pack& gp) {
    io::SingleRead genome("ref", gp.genome.str());
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapRead(genome).path();
    auto path2 = mapper->MapRead(!genome).path();
    return omnigraph::visualization::DefaultColorer(gp.g, path1, path2);
}

inline double l2_norm(const AbundanceVector& v, size_t sample_cnt) {
    double s = 0.;
    for (size_t i = 0; i < sample_cnt; ++i) {
        s += v[i] * v[i];
    }
    return std::sqrt(s);
}

inline double cosine_sim(const AbundanceVector& v1, const AbundanceVector& v2, size_t sample_cnt) {
    double s = 0.;
    for (size_t i = 0; i < sample_cnt; ++i) {
        s += v1[i] * v2[i];
    }
    return s / (l2_norm(v1, sample_cnt) * l2_norm(v2, sample_cnt));
}

template<class Graph>
class EdgeAbundance: public omnigraph::GraphActionHandler<Graph> {
    typedef map<EdgeId, AbundanceVector> Storage;
    typedef Storage::const_iterator const_iterator;
    Storage edge_abundance_;
    const ContigAbundanceCounter& abundance_counter_;

public:
    EdgeAbundance(const Graph& g, const ContigAbundanceCounter& abundance_counter) :
        omnigraph::GraphActionHandler<Graph>(g, "EdgeAbundance"),
        abundance_counter_(abundance_counter){}

    void Fill() {
        for (auto it = this->g().ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            HandleAdd(*it);
        }
    }

    virtual void HandleAdd(EdgeId e) override {
        auto ab = abundance_counter_(this->g().EdgeNucls(e).str());
        if (!ab) {
            INFO("Couldn't estimate abundance of edge " << this->g().str(e));
        } else {
            edge_abundance_[e] = *ab;
        }
    }

    const_iterator begin() const {
        return edge_abundance_.begin();
    }

    const_iterator end() const {
        return edge_abundance_.end();
    }

    const_iterator find(EdgeId e) const {
        return edge_abundance_.find(e);
    }

    size_t count(EdgeId e) const {
        return edge_abundance_.count(e);
    }

private:
    DECL_LOGGER("EdgeAbundance");
};

template<class Graph>
class AggressiveClearing: public omnigraph::EdgeProcessingAlgorithm<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    const size_t sample_cnt_;
    const EdgeAbundance<Graph>& edge_abundance_;
    const AbundanceVector base_profile_;
    const double similarity_threshold_;
    const double norm_ratio_threshold_;
    EdgeRemover<Graph> edge_remover_;
    pred::TypedPredicate<EdgeId> topological_condition_;

protected:
    virtual bool ProcessEdge(EdgeId e) override {
        DEBUG("Processing edge " << this->g().str(e));
        if (!topological_condition_(e)) {
            DEBUG("Topological condition failed");
            return false;
        }
        auto it = edge_abundance_.find(e);
        if (it == edge_abundance_.end()) {
            DEBUG("Edge " << this->g().str(e) << " did not have valid abundance profile");
            return false;
        }
        const auto& profile = it->second;
        DEBUG("Edge profile " << PrintVector(profile, sample_cnt_));
        double sim = cosine_sim(profile, base_profile_, sample_cnt_);
        double norm_ratio = l2_norm(profile, sample_cnt_) / l2_norm(base_profile_, sample_cnt_);

        DEBUG("Similarity between edge and base profiles " << sim);
        DEBUG("Norm ratio " << norm_ratio);
        if (math::ls(norm_ratio, norm_ratio_threshold_) 
                || math::ls(sim, similarity_threshold_)) {
            DEBUG("Removing edge " << this->g().str(e));

            edge_remover_.DeleteEdge(e);
            return true;
        }
        return false;
    }

public:
    AggressiveClearing(Graph &g,
                       size_t sample_cnt,
                       const EdgeAbundance<Graph>& edge_abundance,
                       const AbundanceVector& base_profile,
                       double similarity_threshold,
                       double norm_ratio_threshold,
                       const std::function<void(EdgeId)> &removal_handler = 0) :
        EdgeProcessingAlgorithm<Graph>(g, true),
        sample_cnt_(sample_cnt),
        edge_abundance_(edge_abundance),
        base_profile_(base_profile),
        similarity_threshold_(similarity_threshold),
        norm_ratio_threshold_(norm_ratio_threshold),
        edge_remover_(g, removal_handler),
        topological_condition_(pred::Or(AlternativesPresenceCondition<Graph>(g), TipCondition<Graph>(g))) { 
            DEBUG("Base profile " << PrintVector(base_profile_, sample_cnt_));
        }
private:
    DECL_LOGGER("AggressiveClearing");
};

class SeriesAnalysis : public spades::AssemblyStage {

    boost::optional<AbundanceVector> InferAbundance(const std::string bin_mult_fn,
                                                    const std::string& b_id,
                                                    size_t sample_cnt) const {
        path::CheckFileExistenceFATAL(bin_mult_fn);

        ifstream is(bin_mult_fn);
        vector<AbundanceVector> abundances;
        while (true) {
            string name;
            is >> name;
            if (!is.fail()) {
                AbundanceVector vec;
                for (size_t i = 0; i < sample_cnt; ++i) {
                    is >> vec[i];
                    VERIFY(!is.fail());
                }
                if (name == b_id) {
                    abundances.push_back(vec);
                }
            } else {
                INFO("Read " << abundances.size() << " profiles for bin " << b_id);
                break;
            }
        }

        return boost::optional<AbundanceVector>(MeanVector(abundances, sample_cnt));
    }

public:
    SeriesAnalysis() : AssemblyStage("Series Analysis", "series_analysis") { }

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *) {
        std::string cfg = cfg::get().series_analysis;//"/Sid/snurk/mts/out/infant_gut_2/reassembly.yaml";
        if (cfg.empty()) {
            INFO("No series analysis config was provided via --series-analysis");
            return;
        } else {
            INFO("Series analysis enabled with config " << cfg);
        }
        auto Buf = llvm::MemoryBuffer::getFile(cfg);
        if (!Buf)
            throw std::string("Failed to load config file " + cfg);

        llvm::yaml::Input yin(*Buf.get());
        SeriesAnalysisConfig config;
        yin >> config;
        
        ContigAbundanceCounter abundance_counter(config.k, config.sample_cnt, 
                                                 SingleClusterAnalyzer(config.sample_cnt, 2., 0.4),
                                                 cfg::get().tmp_dir);

        abundance_counter.Init(config.kmer_mult);
        boost::optional<AbundanceVector> bin_profile = InferAbundance(config.bin_prof, config.bin, config.sample_cnt);
        if (!bin_profile) {
            ERROR("Couldn't estimate profile of bin");
            return;
        }

        EdgeAbundance<Graph> edge_abundance(gp.g, abundance_counter);
        edge_abundance.Fill();

        gp.EnsureBasicMapping();
        gp.FillQuality();
        omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
        auto colorer = DefaultGPColorer(gp);
        path::make_dir(cfg::get().output_dir + "pictures/");
        QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
                                       cfg::get().output_dir + "pictures/");

        INFO("Launching aggressive graph clearing");
        //positive quality edges removed (folder colored_edges_deleted)
        AggressiveClearing<Graph> clearing(gp.g, config.sample_cnt, edge_abundance,
                                            *bin_profile, 0.8, 0.3, [&](EdgeId e) {
                        qual_removal_handler.HandleDelete(e);});
        clearing.Run();
        INFO("Graph clearing finished");

        INFO("Drawing edges with failed abundance estimate")
        path::make_dir(cfg::get().output_dir + "pictures_no_ab/");
        QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler2(gp.g, gp.edge_qual, labeler, colorer,
                                       cfg::get().output_dir + "pictures_no_ab/");

        for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            if (edge_abundance.count(e) == 0) {
                qual_removal_handler2.HandleDelete(e);
            } 
        }
    }
};

}

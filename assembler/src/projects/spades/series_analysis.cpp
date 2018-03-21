//***************************************************************************
//* Copyright (c) 2016-2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/handlers/id_track_handler.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "modules/simplification/tip_clipper.hpp"
#include "projects/mts/contig_abundance.hpp"
#include "io/reads/osequencestream.hpp"
#include "series_analysis.hpp"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

namespace debruijn_graph {

struct SeriesAnalysisConfig {
    uint k;
    uint sample_cnt;
    uint frag_size;
    uint min_len;

    std::string kmer_mult, bin, bin_prof, edges_sqn, edges_mpl, edge_fragments_mpl;
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
        io.mapRequired("min_len", cfg.min_len);
        io.mapRequired("edges_sqn", cfg.edges_sqn);
        io.mapRequired("edges_mpl", cfg.edges_mpl);
        io.mapRequired("edge_fragments_mpl", cfg.edge_fragments_mpl);
        io.mapRequired("frag_size", cfg.frag_size);
    }
};

} }

namespace debruijn_graph {

typedef Profile<Abundance> AbundanceVector;
typedef ProfileCounter<Abundance> ContigAbundanceCounter;

template<class graph_pack>
shared_ptr<visualization::graph_colorer::GraphColorer<typename graph_pack::graph_t>> DefaultGPColorer(
    const graph_pack& gp) {
    io::SingleRead genome("ref", gp.genome.str());
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapRead(genome).path();
    auto path2 = mapper->MapRead(!genome).path();
    return visualization::graph_colorer::DefaultColorer(gp.g, path1, path2);
}

inline double l2_norm(const AbundanceVector& v) {
    double s = 0.;
    for (auto val : v) {
        s += val * val;
    }
    return std::sqrt(s);
}

inline double cosine_sim(const AbundanceVector& v1, const AbundanceVector& v2) {
    double s = 0.;
    for (size_t i = 0; i < v1.size(); ++i) {
        s += v1[i] * v2[i];
    }
    return s / (l2_norm(v1) * l2_norm(v2));
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
    const EdgeAbundance<Graph>& edge_abundance_;
    const AbundanceVector base_profile_;
    const double similarity_threshold_;
    const double norm_ratio_threshold_;
    EdgeRemover<Graph> edge_remover_;
    func::TypedPredicate<EdgeId> topological_condition_;

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
        DEBUG("Edge profile " << PrintVector(profile));
        double sim = cosine_sim(profile, base_profile_);
        double norm_ratio = l2_norm(profile) / l2_norm(base_profile_);

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
                       const EdgeAbundance<Graph>& edge_abundance,
                       const AbundanceVector& base_profile,
                       double similarity_threshold,
                       double norm_ratio_threshold,
                       const std::function<void(EdgeId)> &removal_handler = 0) :
        EdgeProcessingAlgorithm<Graph>(g, true),
        edge_abundance_(edge_abundance),
        base_profile_(base_profile),
        similarity_threshold_(similarity_threshold),
        norm_ratio_threshold_(norm_ratio_threshold),
        edge_remover_(g, removal_handler),
        topological_condition_(func::Or(AlternativesPresenceCondition<Graph>(g), TipCondition<Graph>(g))) {
            DEBUG("Base profile " << PrintVector(base_profile_));
        }
private:
    DECL_LOGGER("AggressiveClearing");
};

boost::optional<AbundanceVector> InferAbundance(const std::string& bin_mult_fn,
                                                const std::string& b_id) {
    fs::CheckFileExistenceFATAL(bin_mult_fn);

    ifstream is(bin_mult_fn);
    std::vector<AbundanceVector> abundances;
    std::string name;
    while (true) {
        is >> name;
        if (!is.fail()) {
            if (name != b_id) {
                is.ignore(numeric_limits<std::streamsize>::max(), '\n');
                continue;
            }
            AbundanceVector vec(KmerProfileIndex::SampleCount(), 0.0);
            for (size_t i = 0; i < vec.size(); ++i) {
                is >> vec[i];
                VERIFY(!is.fail());
            }
            abundances.push_back(vec);
        } else {
            INFO("Read " << abundances.size() << " profiles for bin " << b_id);
            break;
        }
    }
    return boost::optional<AbundanceVector>(MeanVector(abundances));
}

void PrintEdgeFragmentProfiles(const conj_graph_pack &gp, const ContigAbundanceCounter &abundance_counter,
                               size_t split_length, size_t min_len, std::ostream &os) {
    for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        io::SingleRead full_contig(std::to_string(gp.g.int_id(e)), gp.g.EdgeNucls(e).str());
        for (size_t i = 0; i < full_contig.size(); i += split_length) {
            if (full_contig.size() - i < min_len) {
                DEBUG("Fragment shorter than min_length_bound " << min_len);
                break;
            }

            io::SingleRead contig = full_contig.Substr(i, std::min(i + split_length, full_contig.size()));

            DEBUG("Processing fragment # " << (i / split_length) << " with id " << contig.name());

            auto abundance_vec = abundance_counter(contig.GetSequenceString(), contig.name());

            if (abundance_vec) {
                size_t len = contig.GetSequenceString().size();
                os << contig.name() << " " << len << " " << PrintVector(*abundance_vec) << std::endl;
                //copy(abundance_vec->begin(), abundance_vec->begin() + config.sample_cnt,
                //     ostream_iterator<Mpl>(ss, " "));
                DEBUG("Successfully estimated abundance of " << contig.name());
            } else {
                DEBUG("Failed to estimate abundance of " << contig.name());
            }
        }
    }
}

void SeriesAnalysis::run(conj_graph_pack &gp, const char *) {
    std::string cfg = cfg::get().series_analysis;
    INFO("Series analysis enabled with config " << cfg);

    auto buf = llvm::MemoryBuffer::getFile(cfg);
    VERIFY_MSG(buf, "Failed to load config file " + cfg);

    llvm::yaml::Input yin(*buf.get());
    SeriesAnalysisConfig config;
    yin >> config;

    KmerProfileIndex::SetSampleCount(config.sample_cnt);

    DEBUG("Initiating abundance counter");
    ContigAbundanceCounter abundance_counter =
        MakeTrivial<Abundance>(config.k, config.kmer_mult);

    DEBUG("Abundance counter ready");

    if (!config.edges_sqn.empty()) {
        io::OFastaReadStream oss(config.edges_sqn);
        for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            string s = gp.g.EdgeNucls(e).str();
            oss << io::SingleRead(io::MakeContigId(gp.g.int_id(e), s.size()), s);
        }
    }

    if (!config.edges_mpl.empty()) {
        ofstream os(config.edges_mpl);
        PrintEdgeFragmentProfiles(gp, abundance_counter, -1ul, config.min_len, os);
    }

    if (!config.edge_fragments_mpl.empty()) {
        ofstream os(config.edge_fragments_mpl);
        PrintEdgeFragmentProfiles(gp, abundance_counter, config.frag_size, config.min_len, os);
    }

//    boost::optional<AbundanceVector> bin_profile = InferAbundance(config.bin_prof, config.bin);
//    if (!bin_profile) {
//        ERROR("Couldn't estimate profile of bin");
//        return;
//    }
//
//    EdgeAbundance<Graph> edge_abundance(gp.g, abundance_counter);
//    edge_abundance.Fill();
//
//    gp.EnsureBasicMapping();
//    gp.FillQuality();
//    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
//    auto colorer = DefaultGPColorer(gp);
//
//    /*
//    fs::make_dir(cfg::get().output_dir + "pictures/");
//    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
//                                   cfg::get().output_dir + "pictures/");
//
//    INFO("Launching aggressive graph clearing");
//    //positive quality edges removed (folder colored_edges_deleted)
//    AggressiveClearing<Graph> clearing(gp.g, edge_abundance,
//                                        *bin_profile, 0.8, 0.3, [&](EdgeId e) {
//                    qual_removal_handler.HandleDelete(e);});
//    clearing.Run();
//    INFO("Graph clearing finished");
//    */
//
//    INFO("Drawing edges with failed abundance estimate")
//    fs::make_dir(cfg::get().output_dir + "pictures_no_ab/");
//    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler2(gp.g, gp.edge_qual, labeler, colorer,
//                                   cfg::get().output_dir + "pictures_no_ab/");
//
//    for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
//        EdgeId e = *it;
//        if (edge_abundance.count(e) == 0) {
//            qual_removal_handler2.HandleDelete(e);
//        }
//    }
}

}

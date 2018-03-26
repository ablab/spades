#include "cloud_scaffold_graph_gap_closer.hpp"
#include "pipeline/extenders_logic.hpp"

namespace path_extend {
shared_ptr<ScaffoldGraphGapCloser> ReadCloudScaffoldGraphGapCloserConstructor::ConstructGraphBasedGapCloser(
    const ReadCloudScaffoldGraphGapCloserConstructor::ScaffoldGraph &graph, size_t edge_length_threshold) const {
    auto path_extender = ConstructExtender(edge_length_threshold);
    auto length_predicate = [this, edge_length_threshold](const EdgeId &edge) {
      return this->gp_.g.length(edge) >= edge_length_threshold;
    };
    ScaffoldGraphExtractor extractor;
    DEBUG("Getting edge map");
    auto long_edge_to_vertex = extractor.GetFirstEdgeMap(graph, length_predicate);
    DEBUG("Got edge map");
    auto tip_extender = make_shared<TipExtender>(gp_.g, path_extender, edge_length_threshold);
    auto tip_searcher =
        make_shared<PathExtenderTipSearcher>(gp_.g, tip_extender, long_edge_to_vertex, edge_length_threshold);
    auto gap_closer = make_shared<TipFinderGapCloser>(tip_searcher);
    return gap_closer;
}
shared_ptr<PathExtender> ReadCloudScaffoldGraphGapCloserConstructor::ConstructExtender(size_t seed_edge_length) const {
    path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                  cfg::get().pe_params,
                                                  cfg::get().ss,
                                                  cfg::get().output_dir,
                                                  cfg::get().mode,
                                                  cfg::get().uneven_depth,
                                                  cfg::get().avoid_rc_connections,
                                                  cfg::get().use_scaffolder);
    auto dataset_info = cfg::get().ds;
    PELaunchSupport support(dataset_info, params);

//    ScaffoldingUniqueEdgeAnalyzer analyzer(gp_, edge_length_threshold, 50.0);
//    analyzer.FillUniqueEdgeStorage(unique_storage);
    //fixme temporary fix to deal with reference fields in LoopDetectingExtender
    ScaffoldingUniqueEdgeStorage *unique_storage = new ScaffoldingUniqueEdgeStorage;
    const GraphCoverageMap *cover_map = new GraphCoverageMap(gp_.g);
    UsedUniqueStorage *used_unique_storage = new UsedUniqueStorage(*unique_storage);

    UniqueData unique_data;
    ExtendersGenerator generator(dataset_info, params, gp_, *cover_map,
                                 unique_data, *used_unique_storage, support);
//    auto basic_extenders = generator.MakeBasicExtenders();
    optional<std::size_t> paired_lib_index;
    for (std::size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info.reads[lib_index];
        if (lib.type() == io::LibraryType::Clouds10x) {
            paired_lib_index = lib_index;
            break;
        }
    }
    VERIFY(paired_lib_index.is_initialized());
    const auto &lib = dataset_info.reads[paired_lib_index.get()];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[paired_lib_index.get()]);
    shared_ptr<WeightCounter> weight_counter = make_shared<ReadCountWeightCounter>(gp_.g, paired_lib);
    auto opts = params.pset.extension_options;
    double weight_threshold = opts.weight_threshold;

    auto pe_extension_chooser = generator.MakeSimpleExtensionChooser(paired_lib_index.get());

    const size_t reliable_edge_length = params_.reliable_edge_length_;
    const size_t tail_threshold = params_.tail_threshold_;
    const size_t distance_bound = params_.distance_bound_;
    const double score_threshold = params_.extender_score_threshold_;
    const double relative_coverage_threshold = params_.relative_coverage_threshold_;

    auto barcode_extractor =
        std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    auto entry_collector = std::make_shared<SimpleBarcodeEntryCollector>(gp_.g, barcode_extractor, reliable_edge_length,
                                                                         seed_edge_length, tail_threshold,
                                                                         relative_coverage_threshold);
    auto read_cloud_extension_chooser = make_shared<ReadCloudExtensionChooser>(gp_.g,
                                                                               weight_counter,
                                                                               weight_threshold,
                                                                               barcode_extractor,
                                                                               entry_collector,
                                                                               score_threshold);
    auto composite_chooser =
        make_shared<CompositeExtensionChooser>(gp_.g, pe_extension_chooser, read_cloud_extension_chooser);

    size_t insert_size = paired_lib->GetISMax();
    bool investigate_short_loops = false;
    bool use_short_loops_cov_resolver = false;

    INFO("Unique check enabled: " << used_unique_storage->UniqueCheckEnabled());
    auto read_cloud_extender = make_shared<ReadCloudExtender>(gp_,
                                                              *cover_map,
                                                              *used_unique_storage,
                                                              composite_chooser,
                                                              insert_size,
                                                              investigate_short_loops,
                                                              use_short_loops_cov_resolver,
                                                              weight_threshold,
                                                              reliable_edge_length,
                                                              tail_threshold,
                                                              distance_bound);
    return read_cloud_extender;
}
ReadCloudScaffoldGraphGapCloserConstructor::ReadCloudScaffoldGraphGapCloserConstructor(const conj_graph_pack &gp_,
                                                                                       const ScaffoldGraphGapCloserParams &params)
    : gp_(gp_), params_(params) {}
shared_ptr<ScaffoldGraphGapCloser> ReadCloudScaffoldGraphGapCloserConstructor::ConstructCloudBasedGapCloser(
        const ReadCloudScaffoldGraphGapCloserConstructor::ScaffoldGraph &graph, size_t edge_length_threshold) const {
    auto path_extender = ConstructExtender(edge_length_threshold);
    auto length_predicate = [this, edge_length_threshold](const EdgeId &edge) {
      return this->gp_.g.length(edge) >= edge_length_threshold;
    };
    auto barcode_extractor =
        std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    auto entry_collector = std::make_shared<SimpleBarcodeEntryCollector>(gp_.g, barcode_extractor,
                                                                         params_.reliable_edge_length_,
                                                                         edge_length_threshold,
                                                                         params_.tail_threshold_,
                                                                         params_.relative_coverage_threshold_);
    ScaffoldGraphExtractor extractor;
    DEBUG("Getting edge map");
    auto long_edge_to_vertex = extractor.GetFirstEdgeMap(graph, length_predicate);
    DEBUG("Got edge map");
    auto tip_extender = make_shared<TipExtender>(gp_.g, path_extender, edge_length_threshold);
    auto gap_closer = make_shared<ScoreFunctionGapCloser>(gp_.g, tip_extender, entry_collector, params_.tip_score_threshold_);
    return gap_closer;
}

ScaffoldGraphGapCloserParams::ScaffoldGraphGapCloserParams(const size_t reliable_edge_length_,
                                                           const size_t tail_threshold_,
                                                           const size_t distance_bound_,
                                                           const double extender_score_threshold_,
                                                           const double tip_score_threshold_,
                                                           const double relative_coverage_threshold_)
    : reliable_edge_length_(reliable_edge_length_),
      tail_threshold_(tail_threshold_),
      distance_bound_(distance_bound_),
      extender_score_threshold_(extender_score_threshold_),
      tip_score_threshold_(tip_score_threshold_),
      relative_coverage_threshold_(relative_coverage_threshold_) {}
}


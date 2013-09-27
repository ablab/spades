//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * distance_estimation.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "late_pair_info_count.hpp"
#include "gap_closer.hpp"
#include "check_tools.hpp"
#include "pair_info_improver.hpp"
#include "long_read_storage.hpp"

#include "de/paired_info.hpp"
#include "de/weighted_distance_estimation.hpp"
#include "de/extensive_distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"

#include <set>

namespace debruijn_graph {

void estimate_with_estimator(const Graph& /*graph*/,
                             const AbstractDistanceEstimator<Graph>& estimator,
                             const PairInfoWeightFilter<Graph>& filter,
                             PairedIndexT& clustered_index)
{
  using debruijn_graph::estimation_mode;
  DEBUG("Estimating distances");

  if (cfg::get().use_multithreading) {
    estimator.EstimateParallel(clustered_index, cfg::get().max_threads);
  } else {
    estimator.Estimate(clustered_index);
  }

  INFO("Filtering info");
  filter.Filter(clustered_index);
  DEBUG("Info Filtered");
}

bool prepare_scaffold_index(conj_graph_pack& gp,
                               const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                               const PairedIndexT& paired_index,
                               PairedIndexT& clustered_index) {

    double is_var = lib.data().insert_size_deviation;
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(
            cfg::get().de.linkage_distance_coeff * is_var);
    GraphDistanceFinder<Graph> dist_finder(gp.g, (size_t) math::round(lib.data().mean_insert_size),
                                           lib.data().read_length, delta);
    size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
    boost::function<double(int)> weight_function;

    DEBUG("Retaining insert size distribution for it");
    if (lib.data().insert_size_distribution.size() == 0) {
        return false;
    }

    WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
    DEBUG("Weight Wrapper Done");
    weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);

    PairInfoWeightFilter<Graph> filter(gp.g, 0.);
    DEBUG("Weight Filter Done");

    const AbstractDistanceEstimator<Graph>& estimator =
            SmoothingDistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
                    weight_function, linkage_distance, max_distance,
                    cfg::get().ade.threshold, cfg::get().ade.range_coeff,
                    cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
                    cfg::get().ade.min_peak_points, cfg::get().ade.inv_density,
                    cfg::get().ade.percentage,
                    cfg::get().ade.derivative_threshold, true);
    estimate_with_estimator(gp.g, estimator, filter, clustered_index);

    return true;
}

void estimate_distance(conj_graph_pack& gp,
                       const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                       const PairedIndexT& paired_index,
                       PairedIndexT& clustered_index,
                       PairedIndexT& scaffold_index)
{
  using debruijn_graph::estimation_mode;

  if (!cfg::get().developer_mode) {
    clustered_index.Attach();
    clustered_index.Init();
  }

  const debruijn_config& config = cfg::get();
  if (config.paired_mode) {

      size_t delta = size_t(lib.data().insert_size_deviation);
    size_t linkage_distance = size_t(config.de.linkage_distance_coeff * lib.data().insert_size_deviation);
    GraphDistanceFinder<Graph> dist_finder(gp.g,  (size_t)math::round(lib.data().mean_insert_size), lib.data().read_length, delta);
    size_t max_distance = size_t(config.de.max_distance_coeff * lib.data().insert_size_deviation);

    boost::function<double(int)> weight_function;

    if (config.est_mode == em_weighted                                        // in these cases we need a weight function
     || config.est_mode == em_smoothing                                       // to estimate graph distances in the
     || config.est_mode == em_extensive)                                      // histogram
    {
      if (lib.data().insert_size_distribution.size() == 0) {
          WARN("No insert size distribution found, stopping distance estimation");
          return;
      }
      WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
      DEBUG("Weight Wrapper Done");
      weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);
    }
    else
      weight_function = UnityFunction;

    PairInfoWeightFilter<Graph> filter(gp.g, config.de.filter_threshold);
    DEBUG("Weight Filter Done");

    switch (config.est_mode)
    {
      case em_simple :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  DistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
                                              linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, filter, clustered_index);
          break;
        }
      case em_weighted :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  WeightedDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, filter, clustered_index);
          break;
        }
      case em_extensive :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  ExtensiveDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, filter, clustered_index);
          break;
        }
      case em_smoothing :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  SmoothingDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance,
                      config.ade.threshold,
                      config.ade.range_coeff,
                      config.ade.delta_coeff, config.ade.cutoff,
                      config.ade.min_peak_points,
                      config.ade.inv_density,
                      config.ade.percentage,
                      config.ade.derivative_threshold);

          estimate_with_estimator(gp.g, estimator, filter, clustered_index);
          break;
        }
    }

    INFO("Refining clustered pair information ");                              // this procedure checks, whether index
    RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
    DEBUG("The refining of clustered pair information has been finished ");    // if so, it resolves such conflicts.

    INFO("Filling paired information");
    PairInfoImprover<Graph> improver(gp.g, clustered_index, lib);
    improver.ImprovePairedInfo(config.use_multithreading, config.max_threads);

    //DE for scaffolds
    INFO("Estimating distances for scaffolding");
    if (!prepare_scaffold_index(gp, lib, paired_index, scaffold_index)) {
        WARN("This lib can not be used for scaffolding");
    }
  }
}

void load_distance_estimation(conj_graph_pack& gp,
                              PairedIndicesT& paired_indices,
                              PairedIndicesT& clustered_indices,
                              PairedIndicesT& scaffold_indices,
                              path::files_t* used_files,
                              PathStorage<Graph>& long_reads,
                              LongReadContainerT& single_long_reads) {
  string p;
  if (cfg::get().entry_point == ws_repeats_resolving && cfg::get().pacbio_test_on)
      p = path::append_path(cfg::get().load_from, "pacbio_aligning");
  else
      p = path::append_path(cfg::get().load_from, "distance_estimation");
  used_files->push_back(p);
  ScanAll(p, gp, paired_indices, clustered_indices, scaffold_indices, single_long_reads);
  if (cfg::get().entry_point == ws_repeats_resolving && cfg::get().pacbio_test_on) {
      INFO(" long reads loading form " <<p);
      long_reads.LoadFromFile(p + ".mpr");
      INFO(" long reads loaded " );
  }
  //TODO: load single long reads
  load_lib_data(p);
}
//
//void load_pacbio_aligned(PathStorage<Graph>& long_reads, path::files_t* used_files)
//{
//    string p = path::append_path(cfg::get().load_from, "pacbio_aligning");
//    used_files->push_back(p);
//    ScanAll(p, gp_, paired_indices, clustered_indices);
//    long_reads.LoadFromFile(p + ".mpr");
//    load_lib_data(p);
//}

void save_distance_estimation(const conj_graph_pack& gp,
                              const PairedIndicesT& paired_indices,
                              const PairedIndicesT& clustered_indices,
                              const PairedIndicesT& scaffold_indices,
                              const LongReadContainerT& single_long_reads)
{
  if (cfg::get().make_saves || (cfg::get().paired_mode && cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles)) {
    string p = path::append_path(cfg::get().output_saves, "distance_estimation");
    INFO("Saving current state to " << p);
    PrintAll(p, gp, paired_indices, clustered_indices, scaffold_indices, single_long_reads);
    write_lib_data(p);
  }
  //TODO: load single long read
}


//void count_estimated_info_stats(const conj_graph_pack& gp,
//                                const PairedIndexT& paired_index,
//                                const PairedIndexT& clustered_index)
//{
//  CountClusteredPairedInfoStats(gp, paired_index, clustered_index);
//}

void exec_distance_estimation(conj_graph_pack& gp,
                              PairedIndicesT& paired_indices,
                              PairedIndicesT& clustered_indices,
                              PairedIndicesT& scaffold_indices,
                              PathStorage<Graph>& pacbio_reads,
                              LongReadContainerT& single_long_reads) {
  if (cfg::get().entry_point <= ws_distance_estimation) {
    exec_late_pair_info_count(gp, paired_indices, single_long_reads);
    if (cfg::get().paired_mode) {
        INFO("STAGE == Estimating Distance");
		for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
	        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0 &&
	                (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||  cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs)) {
	            estimate_distance(gp, cfg::get_writable().ds.reads[i], paired_indices[i], clustered_indices[i], scaffold_indices[i]);
	        }
		}
    }
    save_distance_estimation(gp, paired_indices, clustered_indices, scaffold_indices, single_long_reads);
  }
  else {
    INFO("Loading Distance Estimation");
    path::files_t used_files;
    load_distance_estimation(gp, paired_indices, clustered_indices, scaffold_indices, &used_files, pacbio_reads, single_long_reads);
    link_files_by_prefix(used_files, cfg::get().output_saves);
  }
  if (cfg::get().paired_mode && cfg::get().paired_info_statistics)
      WARN("Paired statistics is currently unavalable");
    //count_estimated_info_stats(gp, paired_index, clustered_index);            // counting stats for paired index (etalon, false positives, etc.)
}



}

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

#include "de/paired_info.hpp"
#include "de/weighted_distance_estimation.hpp"
#include "de/extensive_distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"

#include <set>

namespace debruijn_graph {

typedef set<Point> Histogram;

void estimate_with_estimator(const Graph& graph,
                             const AbstractDistanceEstimator<Graph>& estimator,
                             const PairedInfoNormalizer<Graph>& normalizer,
                             const PairInfoWeightFilter<Graph>& filter,
                             PairedIndexT& clustered_index)
{
  using debruijn_graph::estimation_mode;
  INFO("Estimating distances");
  PairedIndexT raw_clustered_index(graph);

  if (cfg::get().use_multithreading) {
    estimator.EstimateParallel(raw_clustered_index, cfg::get().max_threads);
  } else {
    estimator.Estimate(raw_clustered_index);
  }

  INFO("Normalizing Weights");
  PairedIndexT normalized_index(graph);

    // temporary fix for scaffolding (I hope) due to absolute thresholds in path_extend
  if (cfg::get().est_mode == em_weighted
   || cfg::get().est_mode == em_smoothing
   || cfg::get().est_mode == em_extensive)
  {
    //TODO: add to config
    double coeff = (cfg::get().ds.single_cell ? (10. / 80.) : (0.2 / 3.00) );
    normalizer.FillNormalizedIndex(raw_clustered_index, normalized_index, coeff);
  }
  else
    normalizer.FillNormalizedIndex(raw_clustered_index, normalized_index);

  INFO("Filtering info");
  filter.Filter(normalized_index, clustered_index);
  DEBUG("Info Filtered");
}

void estimate_distance(conj_graph_pack& gp,
                       const PairedIndexT& paired_index,
                             PairedIndexT& clustered_index)
{
  using debruijn_graph::estimation_mode;
  INFO("Here");
  if (!cfg::get().developer_mode) {
    if (!clustered_index.IsAttached()){
    	clustered_index.Attach();
    	clustered_index.Init();
    }
  }
  INFO("after");

  const debruijn_config& config = cfg::get();
  if (config.paired_mode) {
    INFO("STAGE == Estimating Distance");
    double is_var = config.ds.is_var();
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(config.de.linkage_distance_coeff * is_var);
    GraphDistanceFinder<Graph> dist_finder(gp.g, config.ds.IS(), config.ds.RL(), delta);
    size_t max_distance = size_t(config.de.max_distance_coeff * is_var);
    boost::function<double(int)> weight_function;

    if (config.est_mode == em_weighted                                        // in these cases we need a weight function
     || config.est_mode == em_smoothing                                       // to estimate graph distances in the
     || config.est_mode == em_extensive)                                      // histogram
    {
      INFO("Retaining insert size distribution for it");
      std::map<int, size_t> insert_size_hist = config.ds.hist();
      if (insert_size_hist.size() == 0) {
        auto streams = paired_binary_readers(false, 0);
        GetInsertSizeHistogram(streams, gp, config.ds.IS(), config.ds.is_var(),
                               insert_size_hist);
      }
      WeightDEWrapper wrapper(insert_size_hist, config.ds.IS());
      INFO("Weight Wrapper Done");
      weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);
    }
    else
      weight_function = UnityFunction;

    PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
    if (config.ds.single_cell) {                                              // paired info normalization
      normalizing_f = &TrivialWeightNormalization<Graph>;                     // only in the single-cell case,
    } else {                                                                  // in the case of ``multi-cell''
      // todo reduce number of constructor params                             // we use a trivial weight (equal to 1.)
      PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
                                                          config.ds.IS(), config.ds.is_var(), config.ds.RL(),
                                                          gp.k_value, config.ds.avg_coverage());
      normalizing_f = boost::bind(&PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
                                  weight_normalizer, _1, _2, _3);
    }
    PairedInfoNormalizer<Graph> normalizer(normalizing_f);
    INFO("Normalizer Done");

    PairInfoWeightFilter<Graph> filter(gp.g, config.de.filter_threshold);
    INFO("Weight Filter Done");

    switch (config.est_mode)
    {
      case em_simple :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  DistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
                                              linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
      case em_weighted :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  WeightedDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
      case em_extensive :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  ExtensiveDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
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

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
    }

    INFO("Refining clustered pair information");                              // this procedure checks, whether index
    RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
    DEBUG("The refining of clustered pair information has been finished");    // if so, it resolves such conflicts.

    INFO("Filling paired information");
    PairInfoImprover<Graph> improver(gp.g, clustered_index);
    improver.ImprovePairedInfo(config.use_multithreading, config.max_threads);
    //save_distance_filling(gp, paired_index, clustered_index);

    if (cfg::get().developer_mode && cfg::get().pos.late_threading) {
		FillPos(gp, gp.genome, "10");
		FillPos(gp, !gp.genome, "11");
		if (!cfg::get().pos.contigs_for_threading.empty()
			&& FileExists(cfg::get().pos.contigs_for_threading)) {
		  FillPosWithRC(gp, cfg::get().pos.contigs_for_threading, "thr_");
		}

		if (!cfg::get().pos.contigs_to_analyze.empty()
			&& FileExists(cfg::get().pos.contigs_to_analyze)) {
		  FillPosWithRC(gp, cfg::get().pos.contigs_to_analyze, "anlz_");
		}
	}

  }
}

void load_distance_estimation(conj_graph_pack& gp,
                              PairedIndexT& paired_index,
                              PairedIndexT& clustered_index,
                              path::files_t* used_files)
{
  string p = path::append_path(cfg::get().load_from, "distance_estimation");
  used_files->push_back(p);
  ScanAll(p, gp, paired_index, clustered_index);
  load_estimated_params(p);
}

//bool try_load_distance_filling(conj_graph_pack& gp, PairedIndexT& clustered_index,
    //path::files_t* used_files)
//{
  //WARN("trying to load distance filling");
  //string p = path::append_path(cfg::get().load_from, "distance_filling");

  //FILE* file = fopen((p + ".grp").c_str(), "r");
  //if (file == NULL) {
    //return false;
  //}
  //fclose(file);

  //used_files->push_back(p);

  //clustered_index.Clear();
  //ScannerTraits<conj_graph_pack::graph_t>::Scanner scanner(gp.g,
      //gp.int_ids);
  //ScanClusteredIndex(p, scanner, clustered_index);

  //return true;
//}

//void distance_filling(conj_graph_pack& gp, PairedIndexT& paired_index, PairedIndexT& clustered_index) 
//{
    //path::files_t used_files;
    //if (try_load_distance_filling(gp, clustered_index, &used_files)) {

        //link_files_by_prefix(used_files, cfg::get().output_saves);
        //INFO("Distance filling saves detected and loaded");
    //}
//}

void save_distance_estimation(const conj_graph_pack& gp,
                              const PairedIndexT& paired_index,
                              const PairedIndexT& clustered_index)
{
  if (cfg::get().make_saves || (cfg::get().paired_mode && cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles)) {
    string p = path::append_path(cfg::get().output_saves, "distance_estimation");
    INFO("Saving current state to " << p);
    PrintAll(p, gp, paired_index, clustered_index);
    write_estimated_params(p);
  }
}

//void save_distance_filling(const conj_graph_pack& gp,
                           //const PairedIndexT& paired_index,
                           //const PairedIndexT& clustered_index) {
  //if (cfg::get().make_saves || cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles) {
    //string p = path::append_path(cfg::get().output_saves, "distance_filling");
    //INFO("Saving current state to " << p);
    //PrintAll(p, gp, paired_index, clustered_index);
    //write_estimated_params(p);
  //}
//}

void count_estimated_info_stats(const conj_graph_pack& gp,
                                const PairedIndexT& paired_index,
                                const PairedIndexT& clustered_index)
{
  CountClusteredPairedInfoStats(gp, paired_index, clustered_index);
}

void exec_distance_estimation(conj_graph_pack& gp,
                             PairedIndexT& paired_index,
                             PairedIndexT& clustered_index)
{
  if (cfg::get().entry_point <= ws_distance_estimation) {
    exec_late_pair_info_count(gp, paired_index);
    estimate_distance(gp, paired_index, clustered_index);
    save_distance_estimation(gp, paired_index, clustered_index);
  } else {
    INFO("Loading Distance Estimation");
    path::files_t used_files;
    load_distance_estimation(gp, paired_index, clustered_index, &used_files);
    link_files_by_prefix(used_files, cfg::get().output_saves);
  }
  if (cfg::get().paired_mode && cfg::get().paired_info_statistics)
    count_estimated_info_stats(gp, paired_index, clustered_index);            // counting stats for paired index (etalon, false positives, etc.)
}




// ==== NEW ====


void estimate_distance(conj_graph_pack& gp,
                       io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                       const PairedIndexT& paired_index,
                             PairedIndexT& clustered_index)
{
  using debruijn_graph::estimation_mode;

  if (!cfg::get().developer_mode) {
    clustered_index.Attach();
    clustered_index.Init();
  }

  const debruijn_config& config = cfg::get();
  if (config.paired_mode) {
    INFO("STAGE == Estimating Distance");
    double is_var = config.ds.is_var();
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(config.de.linkage_distance_coeff * is_var);
    GraphDistanceFinder<Graph> dist_finder(gp.g, config.ds.IS(), config.ds.RL(), delta);
    size_t max_distance = size_t(config.de.max_distance_coeff * is_var);
    boost::function<double(int)> weight_function;

    if (config.est_mode == em_weighted                                        // in these cases we need a weight function
     || config.est_mode == em_smoothing                                       // to estimate graph distances in the
     || config.est_mode == em_extensive)                                      // histogram
    {
      INFO("Retaining insert size distribution for it");

      if (lib.data().insert_size_distribution.size() == 0) {
        auto streams = paired_binary_readers(lib, false, 0);
        GetInsertSizeHistogram(streams, gp, lib.data().mean_insert_size, lib.data().insert_size_deviation, lib.data().insert_size_distribution);
      }
      WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
      INFO("Weight Wrapper Done");
      weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);
    }
    else
      weight_function = UnityFunction;

    PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
    if (config.ds.single_cell) {                                              // paired info normalization
      normalizing_f = &TrivialWeightNormalization<Graph>;                     // only in the single-cell case,
    } else {                                                                  // in the case of ``multi-cell''
      // todo reduce number of constructor params                             // we use a trivial weight (equal to 1.)
      PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
              lib.data().mean_insert_size, lib.data().insert_size_deviation, lib.data().read_length,
                                                          gp.k_value, config.ds.avg_coverage());
      normalizing_f = boost::bind(&PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
                                  weight_normalizer, _1, _2, _3);
    }
    PairedInfoNormalizer<Graph> normalizer(normalizing_f);
    INFO("Normalizer Done");

    PairInfoWeightFilter<Graph> filter(gp.g, config.de.filter_threshold);
    INFO("Weight Filter Done");

    switch (config.est_mode)
    {
      case em_simple :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  DistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
                                              linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
      case em_weighted :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  WeightedDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
      case em_extensive :
        {
          const AbstractDistanceEstimator<Graph>&
              estimator =
                  ExtensiveDistanceEstimator<Graph>(gp.g, paired_index,
                      dist_finder, weight_function, linkage_distance, max_distance);

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
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

          estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
          break;
        }
    }

    INFO("Refining clustered pair information");                              // this procedure checks, whether index
    RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
    DEBUG("The refining of clustered pair information has been finished");    // if so, it resolves such conflicts.

    INFO("Filling paired information");
    PairInfoImprover<Graph> improver(gp.g, clustered_index);
    improver.ImprovePairedInfo(config.use_multithreading, config.max_threads);
  }
}

void load_distance_estimation(conj_graph_pack& gp,
                              PairedIndicesT& paired_indices,
                              PairedIndicesT& clustered_indices,
                              path::files_t* used_files)
{
  string p = path::append_path(cfg::get().load_from, "distance_estimation");
  used_files->push_back(p);
  ScanAll(p, gp, paired_indices, clustered_indices);
  load_lib_data(p);
}

void save_distance_estimation(const conj_graph_pack& gp,
                              const PairedIndicesT& paired_indices,
                              const PairedIndicesT& clustered_indices)
{
  if (cfg::get().make_saves || (cfg::get().paired_mode && cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles)) {
    string p = path::append_path(cfg::get().output_saves, "distance_estimation");
    INFO("Saving current state to " << p);
    PrintAll(p, gp, paired_indices, clustered_indices);
    write_lib_data(p);
  }
}


//void count_estimated_info_stats(const conj_graph_pack& gp,
//                                const PairedIndexT& paired_index,
//                                const PairedIndexT& clustered_index)
//{
//  CountClusteredPairedInfoStats(gp, paired_index, clustered_index);
//}

void exec_distance_estimation(conj_graph_pack& gp,
        PairedIndicesT& paired_indices,
        PairedIndicesT& clustered_indices)
{
  if (cfg::get().entry_point <= ws_distance_estimation) {
    exec_late_pair_info_count(gp, paired_indices);

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        estimate_distance(gp, cfg::get_writable().ds.reads[i], paired_indices[i], clustered_indices[i]);
    }
    save_distance_estimation(gp, paired_indices, clustered_indices);
  }
  else {
    INFO("Loading Distance Estimation");
    path::files_t used_files;
    load_distance_estimation(gp, paired_indices, clustered_indices, &used_files);
    link_files_by_prefix(used_files, cfg::get().output_saves);
  }
  if (cfg::get().paired_mode && cfg::get().paired_info_statistics)
      WARN("Paired statistics is currently unavalable");
    //count_estimated_info_stats(gp, paired_index, clustered_index);            // counting stats for paired index (etalon, false positives, etc.)
}



}

#include "config_struct.hpp"

#include "config_common.hpp"
#include "openmp_wrapper.h"

#include "io/reader.hpp"

#include "logger/logger.hpp"

#include <string>
#include <vector>

namespace YAML {
template<>
struct convert<io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> > {
  static Node encode(const io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> &rhs) {
      // First, save the "common" stuff
      Node node = convert<io::SequencingLibraryBase>::encode(rhs);

      // Now, save the remaining stuff
      auto const& data = rhs.data();
      node["read length"] = data.read_length;
      node["mean insert size"] = data.mean_insert_size;
      node["insert size deviation"] = data.insert_size_deviation;
      node["median insert size"] = data.median_insert_size;
      node["insert size mad"] = data.insert_size_mad;
      node["average coverage"] = data.average_coverage;
      node["insert size distribution"] = data.insert_size_distribution;

      return node;
  }

  static bool decode(const Node& node, io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> &rhs) {
      // First, load the "common stuff"
      rhs.load(node);

      // Now load the remaining stuff
      auto& data = rhs.data();
      data.read_length = node["RL"].as<size_t>(0);
      if (data.read_length == 0)
          data.read_length = node["read length"].as<size_t>(0);
      data.mean_insert_size = (double) node["IS"].as<size_t>(0);
      if (data.mean_insert_size == 0.0)
          data.mean_insert_size = node["mean insert size"].as<double>(0);
      data.insert_size_deviation = node["insert size deviation"].as<double>(0.0);
      data.median_insert_size = node["median insert size"].as<double>(0.0);
      data.insert_size_mad = node["insert size mad"].as<double>(0.0);
      data.average_coverage = node["average_coverage"].as<double>(0.0);
      data.insert_size_distribution = node["insert size distribution"].as<decltype(data.insert_size_distribution)>(decltype(data.insert_size_distribution)());

      return true;
  }
};
}

namespace debruijn_graph {
static std::string MakeLaunchTimeDirName() {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer, 80, "%m.%d_%H.%M.%S", timeinfo);
  return std::string(buffer);
}

std::string estimated_param_filename(const std::string& prefix) {
  return prefix + "_est_params.info";
}

void load_lib_data(const std::string& prefix) {
  std::string filename = estimated_param_filename(prefix);

  if (!FileExists(filename)) {
      WARN("Estimates params config " << prefix << " does not exist");
  }

  boost::optional<size_t> lib_count;
  load_param(filename, "lib_count", lib_count);
  if (!lib_count || lib_count != cfg::get().ds.reads.lib_count()) {
      WARN("Estimated params file seems to be incorrect");
      return;
  }

  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
      boost::optional<size_t> sizet_val(0);
      boost::optional<double> double_val(0.);
      load_param(filename, "read_length_" + ToString(i), sizet_val);
      if (sizet_val) {
          cfg::get_writable().ds.reads[i].data().read_length = *sizet_val;
      }
      load_param(filename, "insert_size_" + ToString(i), double_val);
      if (double_val) {
          cfg::get_writable().ds.reads[i].data().mean_insert_size = *double_val;
      }
      load_param(filename, "insert_size_deviation_" + ToString(i), double_val);
      if (double_val) {
          cfg::get_writable().ds.reads[i].data().insert_size_deviation = *double_val;
      }
      load_param(filename, "insert_size_median_" + ToString(i), double_val);
      if (double_val) {
          cfg::get_writable().ds.reads[i].data().median_insert_size = *double_val;
      }
      load_param(filename, "insert_size_mad_" + ToString(i), double_val);
      if (double_val) {
          cfg::get_writable().ds.reads[i].data().insert_size_mad = *double_val;
      }
      load_param(filename, "average_coverage_" + ToString(i), double_val);
      if (double_val) {
          cfg::get_writable().ds.reads[i].data().average_coverage = *double_val;
      }

      load_param_map(filename, "histogram_" + ToString(i), cfg::get_writable().ds.reads[i].data().insert_size_distribution);
  }

}

void write_lib_data(const std::string& prefix) {
  std::string filename = estimated_param_filename(prefix);

  cfg::get().ds.reads.save("foo.txt");

  write_param(filename, "lib_count", cfg::get().ds.reads.lib_count());

  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
      write_param(filename, "read_length_" + ToString(i), cfg::get().ds.reads[i].data().read_length);
      write_param(filename, "insert_size_" + ToString(i), cfg::get().ds.reads[i].data().mean_insert_size);
      write_param(filename, "insert_size_deviation_" + ToString(i), cfg::get().ds.reads[i].data().insert_size_deviation);
      write_param(filename, "insert_size_median_" + ToString(i), cfg::get().ds.reads[i].data().median_insert_size);
      write_param(filename, "insert_size_mad_" + ToString(i), cfg::get().ds.reads[i].data().insert_size_mad);
      write_param(filename, "average_coverage_" + ToString(i), cfg::get().ds.reads[i].data().average_coverage);
      write_param_map(filename, "histogram_" + ToString(i), cfg::get().ds.reads[i].data().insert_size_distribution);
  }
}


void load(debruijn_config::simplification::tip_clipper& tc,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(tc.condition, pt, "condition");
}

void load(working_stage& entry_point,
          boost::property_tree::ptree const& pt, std::string const& key,
          bool complete) {
  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    entry_point = debruijn_config::working_stage_id(ep);
  }
}

void load(resolving_mode& rm, boost::property_tree::ptree const& pt,
          std::string const& key, bool complete) {
  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    rm = debruijn_config::resolving_mode_id(ep);
  }
}

inline void load(construction_mode& con_mode,
        boost::property_tree::ptree const& pt, std::string const& key,
        bool complete) {
    if (complete || pt.find(key) != pt.not_found()) {
        std::string ep = pt.get<std::string>(key);
        con_mode = debruijn_config::construction_mode_id(ep);
    }
}

inline void load(debruijn_config::construction::early_tip_clipper& etc,
        boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;
    load(etc.enable, pt, "enable");
    etc.length_bound = pt.get_optional<size_t>("length_bound");
}

inline void load(debruijn_config::construction& con,
        boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(con.con_mode, pt, "mode", complete);
    load(con.keep_perfect_loops, pt, "keep_perfect_loops", complete);
    load(con.early_tc, pt, "early_tip_clipper", complete);
}

void load(estimation_mode& est_mode,
          boost::property_tree::ptree const& pt, std::string const& key,
          bool complete) {
  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    est_mode = debruijn_config::estimation_mode_id(ep);
  }
}

void load(debruijn_config::simplification::bulge_remover& br,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(br.max_bulge_length_coefficient, pt, "max_bulge_length_coefficient");
  load(br.max_additive_length_coefficient, pt,
       "max_additive_length_coefficient");
  load(br.max_coverage, pt, "max_coverage");
  load(br.max_relative_coverage, pt, "max_relative_coverage");
  load(br.max_delta, pt, "max_delta");
  load(br.max_relative_delta, pt, "max_relative_delta");
}

void load(debruijn_config::simplification::topology_tip_clipper& ttc,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(ttc.length_coeff, pt, "length_coeff");
  load(ttc.plausibility_length, pt, "plausibility_length");
  load(ttc.uniqueness_length, pt, "uniqueness_length");
}

void load(debruijn_config::simplification::relative_coverage_ec_remover& rec,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(rec.max_ec_length_coefficient, pt, "max_ec_length_coefficient");
  load(rec.max_coverage_coeff, pt, "max_coverage_coeff");
  load(rec.coverage_gap, pt, "coverage_gap");
}

void load(debruijn_config::simplification::isolated_edges_remover& ier,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(ier.max_length, pt, "max_length");
  load(ier.max_coverage, pt, "max_coverage");
  load(ier.max_length_any_cov, pt, "max_length_any_cov");
}

void load(debruijn_config::simplification::complex_bulge_remover& cbr,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(cbr.enabled, pt, "enabled");
  load(cbr.pics_enabled, pt, "pics_enabled");
  load(cbr.folder, pt, "folder");
  load(cbr.max_relative_length, pt, "max_relative_length");
  load(cbr.max_length_difference, pt, "max_length_difference");
}

void load(debruijn_config::simplification::erroneous_connections_remover& ec,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(ec.condition, pt, "condition");
}

void load(debruijn_config::simplification::topology_based_ec_remover& tec,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(tec.max_ec_length_coefficient, pt, "max_ec_length_coefficient");
  load(tec.plausibility_length, pt, "plausibility_length");
  load(tec.uniqueness_length, pt, "uniqueness_length");
}

void load(debruijn_config::simplification::interstrand_ec_remover &isec,
          boost::property_tree::ptree const &pt, bool /*complete*/) {
  using config_common::load;
  load(isec.max_ec_length_coefficient, pt, "max_ec_length_coefficient");
  load(isec.uniqueness_length, pt, "uniqueness_length");
  load(isec.span_distance, pt, "span_distance");
}

void load(debruijn_config::simplification::tr_based_ec_remover &trec,
          boost::property_tree::ptree const &pt, bool /*complete*/) {
  using config_common::load;
  load(trec.max_ec_length_coefficient, pt, "max_ec_length_coefficient");
  load(trec.unreliable_coverage, pt, "unreliable_coverage");
  load(trec.uniqueness_length, pt, "uniqueness_length");
}

void load(debruijn_config::simplification::max_flow_ec_remover& mfec,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(mfec.enabled, pt, "enabled");
  load(mfec.max_ec_length_coefficient, pt, "max_ec_length_coefficient");
  load(mfec.plausibility_length, pt, "plausibility_length");
  load(mfec.uniqueness_length, pt, "uniqueness_length");
}

void load(debruijn_config::distance_estimator& de,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(de.linkage_distance_coeff, pt, "linkage_distance_coeff");
  load(de.max_distance_coeff, pt, "max_distance_coeff");
  load(de.filter_threshold, pt, "filter_threshold");
}

void load(debruijn_config::smoothing_distance_estimator& ade,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(ade.threshold, pt, "threshold");
  load(ade.range_coeff, pt, "range_coeff");
  load(ade.delta_coeff, pt, "delta_coeff");
  load(ade.percentage, pt, "percentage");
  load(ade.cutoff, pt, "cutoff");
  load(ade.min_peak_points, pt, "min_peak_points");
  load(ade.inv_density, pt, "inv_density");
  load(ade.derivative_threshold, pt, "derivative_threshold");
}

void load(debruijn_config::repeat_resolver& rr,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(rr.symmetric_resolve, pt, "symmetric_resolve");
  load(rr.mode, pt, "mode");
  load(rr.inresolve_cutoff_proportion, pt, "inresolve_cutoff_proportion");
  load(rr.near_vertex, pt, "near_vertex");
  load(rr.max_distance, pt, "max_distance");
  load(rr.max_repeat_length, pt, "max_repeat_length");
  load(rr.kill_loops, pt, "kill_loops");
}

void load(debruijn_config::coverage_based_rr& cbrr,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

    load(cbrr.coverage_threshold_one_list, pt, "coverage_threshold_one_list");
    load(cbrr.coverage_threshold_match, pt, "coverage_threshold_match");
    load(cbrr.coverage_threshold_global, pt, "coverage_threshold_global");
    load(cbrr.tandem_ratio_lower_threshold, pt, "tandem_ratio_lower_threshold");
    load(cbrr.tandem_ratio_upper_threshold, pt, "tandem_ratio_upper_threshold");
    load(cbrr.repeat_length_upper_threshold, pt, "repeat_length_upper_threshold");

}

void load(debruijn_config::pacbio_processor& pb,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(pb.pacbio_reads, pt, "pacbio_reads");
  load(pb.pacbio_k, pt, "pacbio_k");
  load(pb.additional_debug_info, pt, "additional_debug_info");
  load(pb.pacbio_optimized_sw, pt, "pacbio_optimized_sw");
  load(pb.compression_cutoff, pt, "compression_cutoff");
  load(pb.domination_cutoff, pt, "domination_cutoff");
  load(pb.path_limit_stretching, pt, "path_limit_stretching");
  load(pb.path_limit_pressing, pt, "path_limit_pressing");
  load(pb.gap_closing_iterations, pt, "gap_closing_iterations");
  load(pb.long_seq_limit, pt, "long_seq_limit");
  load(pb.split_cutoff, pt, "split_cutoff");
  load(pb.match_value, pt, "match_value");
  load(pb.mismatch_penalty, pt, "mismatch_penalty");
  load(pb.insertion_penalty, pt, "insertion_penalty");
  load(pb.deletion_penalty, pt, "deletion_penalty");
}


void load(debruijn_config::position_handler& pos,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(pos.max_single_gap, pt, "max_single_gap");
  load(pos.contigs_for_threading, pt, "contigs_for_threading");
  load(pos.contigs_to_analyze, pt, "contigs_to_analyze");
  load(pos.late_threading, pt, "late_threading");
  load(pos.careful_labeling, pt, "careful_labeling");
}

void load(debruijn_config::gap_closer& gc,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(gc.minimal_intersection, pt, "minimal_intersection");
  load(gc.before_simplify, pt, "before_simplify");
  load(gc.in_simplify, pt, "in_simplify");
  load(gc.after_simplify, pt, "after_simplify");
  load(gc.weight_threshold, pt, "weight_threshold");
}

void load(debruijn_config::graph_read_corr_cfg& graph_read_corr,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(graph_read_corr.enable, pt, "enable");
  load(graph_read_corr.output_dir, pt, "output_dir");
  load(graph_read_corr.binary, pt, "binary");
}

void load(debruijn_config::dataset& ds,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(ds.reads_filename, pt, "reads");
  load(ds.single_cell, pt, "single_cell");

  ds.reference_genome_filename = "";
  boost::optional<std::string> refgen =
      pt.get_optional<std::string>("reference_genome");
  if (refgen && *refgen != "N/A") {
    ds.reference_genome_filename = *refgen;
  }
}

void load_reads(debruijn_config::dataset& ds,
        std::string input_dir) {
  if (ds.reads_filename[0] != '/')
    ds.reads_filename = input_dir + ds.reads_filename;
  CheckFileExistenceFATAL(ds.reads_filename);
  ds.reads.load(ds.reads_filename);
}

void load_reference_genome(debruijn_config::dataset& ds,
                           std::string input_dir) {
  if (ds.reference_genome_filename == "") {
    ds.reference_genome = Sequence();
    return;
  }
  if (ds.reference_genome_filename[0] != '/')
    ds.reference_genome_filename = input_dir + ds.reference_genome_filename;
  CheckFileExistenceFATAL(ds.reference_genome_filename);
  io::Reader genome_stream(ds.reference_genome_filename);
  io::SingleRead genome;
  genome_stream >> genome;
  VERIFY(genome.IsValid());
  ds.reference_genome = genome.sequence();
}

void load(debruijn_config::simplification& simp,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(simp.tc, pt, "tc", complete); // tip clipper:
  load(simp.ttc, pt, "ttc", complete); // topology tip clipper:
  load(simp.br, pt, "br", complete); // bulge remover:
  load(simp.ec, pt, "ec", complete); // erroneous connections remover:
  load(simp.rec, pt, "rec", complete); // relative coverage erroneous connections remover:
  load(simp.tec, pt, "tec", complete); // topology aware erroneous connections remover:
  load(simp.trec, pt, "trec", complete); // topology and reliability based erroneous connections remover:
  load(simp.isec, pt, "isec", complete); // interstrand erroneous connections remover (thorn remover):
  load(simp.mfec, pt, "mfec", complete); // max flow erroneous connections remover:
  load(simp.ier, pt, "ier", complete); // isolated edges remover
  load(simp.cbr, pt, "cbr", complete); // complex bulge remover
}

void load(debruijn_config::info_printer& printer,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(printer.print_stats, pt, "print_stats", complete);
  load(printer.write_components, pt, "write_components", complete);
  load(printer.components_for_kmer, pt, "components_for_kmer", complete);
  load(printer.components_for_genome_pos, pt, "components_for_genome_pos",
       complete);
  load(printer.write_components_along_genome, pt,
       "write_components_along_genome", complete);
  load(printer.write_components_along_contigs, pt,
       "write_components_along_contigs", complete);
  load(printer.save_full_graph, pt, "save_full_graph", complete);
  load(printer.write_full_graph, pt, "write_full_graph", complete);
  load(printer.write_full_nc_graph, pt, "write_full_nc_graph", complete);
  load(printer.write_error_loc, pt, "write_error_loc", complete);
}

void load(debruijn_config::info_printers_t& printers,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  using details::info_printer_pos_name;

  debruijn_config::info_printer def;
  load(def, pt, info_printer_pos_name(ipp_default), true);

  for (size_t pos = ipp_default + 1; pos != ipp_total; ++pos) {
    debruijn_config::info_printer printer(def);
    load(printer, pt, info_printer_pos_name(pos), false);

    printers[info_printer_pos(pos)] = printer;
  }
}

// main debruijn config load function
void load(debruijn_config& cfg, boost::property_tree::ptree const& pt,
          bool /*complete*/) {
  using config_common::load;

  load(cfg.K, pt, "K");

  // input options:
  load(cfg.dataset_file, pt, "dataset");
  // input dir is based on dataset file location (all pathes in datasets are relative to its location)
  cfg.input_dir = path::parent_path(cfg.dataset_file);
  if (cfg.input_dir[cfg.input_dir.length() - 1] != '/')
    cfg.input_dir += '/';

  load(cfg.output_base, pt, "output_base");
  if (cfg.output_base[cfg.output_base.length() - 1] != '/')
    cfg.output_base += '/';

  // TODO: remove this option
  load(cfg.run_mode, pt, "run_mode");

  if (cfg.run_mode) {
    load(cfg.project_name, pt, "project_name");
    cfg.output_root =
        cfg.project_name.empty() ?
        (cfg.output_base + "/K" + ToString(cfg.K) + "/") :
        (cfg.output_base + cfg.project_name + "/K"
         + ToString(cfg.K) + "/");
    cfg.output_suffix = MakeLaunchTimeDirName() + "/";
    cfg.output_dir = cfg.output_root + cfg.output_suffix;
  } else {
    cfg.output_root = cfg.output_base + "/K" + ToString(cfg.K) + "/";

    cfg.output_dir = cfg.output_root;
  }

  cfg.output_saves = cfg.output_dir + "saves/";

  load(cfg.log_filename, pt, "log_filename");

  load(cfg.developer_mode, pt, "developer_mode");
  if (cfg.developer_mode) {
    load(cfg.make_saves, pt, "make_saves");
    load(cfg.output_pictures, pt, "output_pictures");
    load(cfg.output_nonfinal_contigs, pt, "output_nonfinal_contigs");
    load(cfg.compute_paths_number, pt, "compute_paths_number");
  } else {
    cfg.make_saves = false;
    cfg.output_pictures = false;
    cfg.output_nonfinal_contigs = false;
    cfg.compute_paths_number = false;
  }
  if (!cfg.make_saves) {
    load(cfg.make_saves, pt, "force_make_saves");
  }

  load(cfg.load_from, pt, "load_from");
  if (cfg.load_from[0] != '/') { // relative path
    if (cfg.run_mode)
      cfg.load_from = cfg.output_root + cfg.load_from;
    else
      cfg.load_from = cfg.output_dir + cfg.load_from;
  }

  load(cfg.entry_point, pt, "entry_point");

  load(cfg.use_additional_contigs, pt, "use_additional_contigs");
  load(cfg.topology_simplif_enabled, pt, "topology_simplif_enabled");
  load(cfg.use_unipaths, pt, "use_unipaths");

  load(cfg.pacbio_test_on, pt, "pacbio_test_on");
  load(cfg.coverage_based_rr_on, pt, "coverage_based_rr_on");
  if (cfg.coverage_based_rr_on) {
    load (cfg.cbrr, pt, "coverage_based_rr");
}
  if (cfg.pacbio_test_on) {
    load(cfg.pb, pt, "pacbio_processor");
  } else {
  }

  load(cfg.additional_contigs, pt, "additional_contigs");

  load(cfg.paired_mode, pt, "paired_mode");
  load(cfg.long_single_mode, pt, "long_single_mode");
  load(cfg.divide_clusters, pt, "divide_clusters");
  load(cfg.mismatch_careful, pt, "mismatch_careful");
  load(cfg.correct_mismatches, pt, "correct_mismatches");
  load(cfg.paired_info_statistics, pt, "paired_info_statistics");
  load(cfg.paired_info_scaffolder, pt, "paired_info_scaffolder");
  load(cfg.cut_bad_connections, pt, "cut_bad_connections");
  load(cfg.componential_resolve, pt, "componential_resolve");
  load(cfg.gap_closer_enable, pt, "gap_closer_enable");

  load(cfg.buffer_size, pt, "buffer_size");
  cfg.buffer_size <<= 20; //turn MB to bytes

  load(cfg.temp_bin_reads_dir, pt, "temp_bin_reads_dir");
  if (cfg.temp_bin_reads_dir[cfg.temp_bin_reads_dir.length() - 1] != '/')
    cfg.temp_bin_reads_dir += '/';
  cfg.temp_bin_reads_path =
      cfg.project_name.empty() ?
      (cfg.output_base + "/" + cfg.temp_bin_reads_dir) :
      (cfg.output_base + cfg.project_name + "/"
       + cfg.temp_bin_reads_dir);

  cfg.temp_bin_reads_info = cfg.temp_bin_reads_path + "INFO";

  cfg.paired_read_prefix = cfg.temp_bin_reads_path + "_paired";
  cfg.single_read_prefix = cfg.temp_bin_reads_path + "_single";

  load(cfg.use_multithreading, pt, "use_multithreading");
  load(cfg.max_threads, pt, "max_threads");
  // Fix number of threads according to OMP capabilities.
  cfg.max_threads = std::min(cfg.max_threads, (size_t) omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads((int) cfg.max_threads);

  load(cfg.max_memory, pt, "max_memory");

  CheckFileExistenceFATAL(cfg.dataset_file);
  boost::property_tree::ptree ds_pt;
  boost::property_tree::read_info(cfg.dataset_file, ds_pt);
  load(cfg.ds, ds_pt, true);

  load(cfg.ade, pt, (cfg.ds.single_cell ? "sc_ade" : "usual_ade")); // advanced distance estimator:
  load(cfg.rr, pt, (cfg.ds.single_cell ? "sc_rr" : "usual_rr")); // repeat resolver:
  load(cfg.pos, pt, "pos"); // position handler:

  load(cfg.est_mode, pt, "estimation_mode");

  load(cfg.rm, pt, "resolving_mode");

  if (cfg.rm == rm_path_extend) {
      load(cfg.de, pt, (cfg.ds.single_cell ? "sc_de" : "usual_de"));
  }
  else {
      load(cfg.de, pt, (cfg.ds.single_cell ? "old_sc_de" : "old_usual_de"));
  }

  cfg.pe_params.name = cfg.ds.single_cell ? "singlecell" : "multicell";
  load(cfg.pe_params, pt, "andrey_params");
  if (!cfg.developer_mode) {
      cfg.pe_params.debug_output = false;
      cfg.pe_params.viz.DisableAll();
      cfg.pe_params.output.DisableAll();
  }
  load(cfg.use_scaffolder, pt, "use_scaffolder");
  if (!cfg.use_scaffolder) {
      cfg.pe_params.param_set.scaffolder_options.on = false;
  }
  load(cfg.avoid_rc_connections, pt, "avoid_rc_connections");

  load(cfg.mask_all, pt, "mask_all");

  load(cfg.con, pt, "construction");
  load(cfg.gc, pt, "gap_closer");
  load(cfg.graph_read_corr, pt, "graph_read_corr");
  load(cfg.need_consensus, pt, "need_consensus");
  load(cfg.uncorrected_reads, pt, "uncorrected_reads");
  load(cfg.mismatch_ratio, pt, "mismatch_ratio");

  load(cfg.simp, pt, "default");

  if (cfg.ds.single_cell)
    load(cfg.simp, pt, "sc", false);

  if (cfg.mismatch_careful)
    load(cfg.simp, pt, "careful", false);

  cfg.simp.cbr.folder = cfg.output_dir + cfg.simp.cbr.folder + "/";

  load(cfg.info_printers, pt, "info_printers");

  load_reads(cfg.ds, cfg.input_dir);
  load_reference_genome(cfg.ds, cfg.input_dir);
}

void load(debruijn_config& cfg, const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  load(cfg, pt, true);
}

};

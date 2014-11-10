#include "config_struct.hpp"

#include "config_common.hpp"
#include "openmp_wrapper.h"

#include "logger/logger.hpp"

#include "io/file_reader.hpp"

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
      node["read length"]                = data.read_length;
      node["average read length"]        = data.avg_read_length;
      node["insert size mean"]           = data.mean_insert_size;
      node["insert size deviation"]      = data.insert_size_deviation;
      node["insert size left quantile"]  = data.insert_size_left_quantile;
      node["insert size right quantile"] = data.insert_size_right_quantile;
      node["insert size median"]         = data.median_insert_size;
      node["insert size mad"]            = data.insert_size_mad;
      node["insert size distribution"]   = data.insert_size_distribution;
      node["average coverage"]           = data.average_coverage;
      node["pi threshold"]               = data.pi_threshold;

      return node;
  }

  static bool decode(const Node& node, io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> &rhs) {
      // First, load the "common stuff"
      rhs.load(node);

      // Now load the remaining stuff
      auto& data = rhs.data();
      data.read_length                = node["read length"].as<size_t>(0);
      data.avg_read_length            = node["average read length"].as<double>(0.0);
      data.mean_insert_size           = node["insert size mean"].as<double>(0.0);
      data.insert_size_deviation      = node["insert size deviation"].as<double>(0.0);
      data.insert_size_left_quantile  = node["insert size left quantile"].as<double>(0.0);
      data.insert_size_right_quantile = node["insert size right quantile"].as<double>(0.0);
      data.median_insert_size         = node["insert size median"].as<double>(0.0);
      data.insert_size_mad            = node["insert size mad"].as<double>(0.0);
      data.insert_size_distribution   = node["insert size distribution"].as<decltype(data.insert_size_distribution)>(decltype(data.insert_size_distribution)());

      data.average_coverage           = node["average coverage"].as<double>(0.0);
      data.pi_threshold               = node["pi threshold"].as<double>(0.0);

      return true;
  }
};

template<>
struct convert<debruijn_graph::debruijn_config::dataset> {
  static Node encode(const debruijn_graph::debruijn_config::dataset &rhs) {
      Node node;

      node["reads"] = rhs.reads;
      node["max read length"] = rhs.RL();
      node["avg read length"] = rhs.aRL();
      node["average coverage"] = rhs.avg_coverage();

      return node;
  }

  static bool decode(const Node& node, debruijn_graph::debruijn_config::dataset &rhs) {
      rhs.set_RL(node["max read length"].as<size_t>(0));
      rhs.set_aRL(node["avg read length"].as<double>(0.0));
      rhs.set_avg_coverage(node["average coverage"].as<double>(0.0));
      rhs.reads = node["reads"];

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

void load_lib_data(const std::string& prefix) {
  // First, load the data into separate libs
  cfg::get_writable().ds.reads.load(prefix + ".lib_data");

  // Now, infer the common parameters
  const auto& reads = cfg::get().ds.reads;
  size_t max_rl = 0;
  double avg_cov = 0.0, avg_rl;
  for (auto it = reads.library_begin(), et = reads.library_end(); it != et; ++it) {
      auto const& data = it->data();
      if (it->is_graph_contructable())
          max_rl = std::max(max_rl, data.read_length);
      if (data.average_coverage > 0)
          avg_cov = data.average_coverage;
      if (data.avg_read_length > 0)
          avg_rl = data.avg_read_length;
  }

  cfg::get_writable().ds.set_RL(max_rl);
  cfg::get_writable().ds.set_aRL(avg_rl);
  cfg::get_writable().ds.set_avg_coverage(avg_cov);
}

void write_lib_data(const std::string& prefix) {
  cfg::get().ds.reads.save(prefix + ".lib_data");
}

void load(debruijn_config::simplification::tip_clipper& tc,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(tc.condition, pt, "condition");
}

void load(resolving_mode& rm, boost::property_tree::ptree const& pt,
          std::string const& key, bool complete) {
  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    rm = debruijn_config::resolving_mode_id(ep);
  }
}

void load(single_read_resolving_mode& rm, boost::property_tree::ptree const& pt,
          std::string const& key, bool complete) {
  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    rm = debruijn_config::single_read_resolving_mode_id(ep);
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
    load(con.read_buffer_size, pt, "read_buffer_size", complete);
    con.read_buffer_size *= 1024 * 1024;
    load(con.early_tc, pt, "early_tip_clipper", complete);
}

inline void load(debruijn_config::sensitive_mapper& sensitive_map,
        boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(sensitive_map.k, pt, "k", complete);
}

inline void load(estimation_mode& est_mode,
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

  load(br.enabled                           , pt,   "enabled"					);
  load(br.max_bulge_length_coefficient		, pt,   "max_bulge_length_coefficient");
  load(br.max_additive_length_coefficient	, pt,
       "max_additive_length_coefficient");
  load(br.max_coverage,                     pt,     "max_coverage");
  load(br.max_relative_coverage,            pt,     "max_relative_coverage");
  load(br.max_delta,                        pt,     "max_delta");
  load(br.max_relative_delta,               pt,     "max_relative_delta");
}

void load(debruijn_config::simplification::topology_tip_clipper& ttc,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(ttc.length_coeff, pt, "length_coeff");
  load(ttc.plausibility_length, pt, "plausibility_length");
  load(ttc.uniqueness_length, pt, "uniqueness_length");
}

void load(debruijn_config::simplification::relative_coverage_comp_remover& rcc,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;
  load(rcc.enabled, pt, "enabled", complete);
  load(rcc.coverage_gap, pt, "coverage_gap", complete);
  load(rcc.length_coeff, pt, "max_length_coeff", complete);
  load(rcc.tip_allowing_length_coeff, pt, "max_length_with_tips_coeff", complete);
  load(rcc.vertex_count_limit, pt, "max_vertex_cnt", complete);
  load(rcc.max_ec_length_coefficient, pt, "max_ec_length_coefficient", complete);
  load(rcc.max_coverage_coeff, pt, "max_coverage_coeff", complete);
}

void load(debruijn_config::simplification::isolated_edges_remover& ier,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(ier.max_length, pt, "max_length");
  load(ier.max_coverage, pt, "max_coverage");
  load(ier.max_length_any_cov, pt, "max_length_any_cov");
}

void load(debruijn_config::simplification::complex_bulge_remover& cbr,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(cbr.enabled, pt, "enabled");
  load(cbr.max_relative_length, pt, "max_relative_length", complete);
  load(cbr.max_length_difference, pt, "max_length_difference", complete);
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

void load(debruijn_config::simplification::hidden_ec_remover& her,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(her.enabled, pt, "enabled");
  load(her.uniqueness_length, pt, "uniqueness_length");
  load(her.unreliability_threshold, pt, "unreliability_threshold");
  load(her.relative_threshold, pt, "relative_threshold");
}

void load(debruijn_config::distance_estimator& de,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(de.linkage_distance_coeff, pt, "linkage_distance_coeff");
  load(de.max_distance_coeff, pt, "max_distance_coeff");
  load(de.max_distance_coeff_scaff, pt, "max_distance_coeff_scaff");
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

inline void load(debruijn_config::ambiguous_distance_estimator& amde,
        boost::property_tree::ptree const& pt, bool){
    using config_common::load;

    load(amde.enabled,						pt,		"enabled");
    load(amde.haplom_threshold,				pt,		"haplom_threshold");
    load(amde.relative_length_threshold,	pt,		"relative_length_threshold");
    load(amde.relative_seq_threshold,		pt,		"relative_seq_threshold");
}

void load(debruijn_config::pacbio_processor& pb,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(pb.pacbio_k, pt, "pacbio_k");
  load(pb.additional_debug_info, pt, "additional_debug_info");
  load(pb.compression_cutoff, pt, "compression_cutoff");
  load(pb.domination_cutoff, pt, "domination_cutoff");
  load(pb.path_limit_stretching, pt, "path_limit_stretching");
  load(pb.path_limit_pressing, pt, "path_limit_pressing");
  load(pb.ignore_middle_alignment, pt, "ignore_middle_alignment");
  load(pb.long_seq_limit, pt, "long_seq_limit");
  load(pb.pacbio_min_gap_quantity, pt, "pacbio_min_gap_quantity");
  load(pb.contigs_min_gap_quantity, pt, "contigs_min_gap_quantity");
}


void load(debruijn_config::position_handler& pos,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(pos.max_mapping_gap, pt, "max_mapping_gap");
  load(pos.max_gap_diff, pt, "max_gap_diff");
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

void load(debruijn_config::kmer_coverage_model& kcm,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(kcm.probability_threshold, pt, "probability_threshold");
  load(kcm.strong_probability_threshold, pt, "strong_probability_threshold");
  load(kcm.coverage_threshold, pt, "coverage_threshold");
  load(kcm.use_coverage_threshold, pt, "use_coverage_threshold");
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
  path::CheckFileExistenceFATAL(ds.reads_filename);
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
  path::CheckFileExistenceFATAL(ds.reference_genome_filename);
  io::FileReadStream genome_stream(ds.reference_genome_filename);
  io::SingleRead genome;
  genome_stream >> genome;
  VERIFY(genome.IsValid());
  ds.reference_genome = genome.sequence();
}

void load(debruijn_config::simplification::presimplification& presimp,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(presimp.enabled, pt, "enabled", complete);
  load(presimp.parallel, pt, "parallel", complete);
  load(presimp.tip_condition, pt, "tip_condition", complete); // pre tip clipper:
  load(presimp.ec_condition, pt, "ec_condition", complete); // pre ec remover:
  load(presimp.ier, pt, "ier", complete);
}

void load(debruijn_config::simplification& simp,
          boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(simp.topology_simplif_enabled, pt, "topology_simplif_enabled", complete);
  load(simp.tc, pt, "tc", complete); // tip clipper:
  load(simp.ttc, pt, "ttc", complete); // topology tip clipper:
  load(simp.br, pt, "br", complete); // bulge remover:
  load(simp.ec, pt, "ec", complete); // erroneous connections remover:
  load(simp.rcc, pt, "rcc", complete); // relative coverage component remover:
  load(simp.tec, pt, "tec", complete); // topology aware erroneous connections remover:
  load(simp.trec, pt, "trec", complete); // topology and reliability based erroneous connections remover:
  load(simp.isec, pt, "isec", complete); // interstrand erroneous connections remover (thorn remover):
  load(simp.mfec, pt, "mfec", complete); // max flow erroneous connections remover:
  load(simp.ier, pt, "ier", complete); // isolated edges remover
  load(simp.cbr, pt, "cbr", complete); // complex bulge remover
  load(simp.her, pt, "her", complete); // hidden ec remover
  load(simp.fast_features, pt, "fast_features", complete); // master switch for speed-up tricks
  load(simp.fast_activation_cov, pt, "fast_activation_cov", complete);
  load(simp.presimp, pt, "presimp", complete); // presimplification
  load(simp.persistent_cycle_iterators, pt, "persistent_cycle_iterators", complete);
  load(simp.disable_br_in_cycle, pt, "disable_br_in_cycle", complete);
//  load(simp.stats_mode, pt, "stats_mode", complete); // temporary stats counting mode
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

void clear(debruijn_config::info_printer& printer) {
    printer.print_stats = false;
    printer.write_components = false;
    printer.components_for_kmer = "";
    printer.components_for_genome_pos = "";
    printer.write_components_along_genome = false;
    printer.save_full_graph = false;
    printer.write_full_graph = false;
    printer.write_full_nc_graph = false;
    printer.write_error_loc = false;
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
  // input dir is based on dataset file location (all paths in datasets are relative to its location)
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
    load(cfg.output_pictures, pt, "output_pictures");
    load(cfg.output_nonfinal_contigs, pt, "output_nonfinal_contigs");
    load(cfg.compute_paths_number, pt, "compute_paths_number");
  } else {
    cfg.output_pictures = false;
    cfg.output_nonfinal_contigs = false;
    cfg.compute_paths_number = false;
  }

  load(cfg.load_from, pt, "load_from");
  if (cfg.load_from[0] != '/') { // relative path
    if (cfg.run_mode)
      cfg.load_from = cfg.output_root + cfg.load_from;
    else
      cfg.load_from = cfg.output_dir + cfg.load_from;
  }

  load(cfg.tmp_dir, pt, "tmp_dir");
  if (cfg.tmp_dir[0] != '/') { // relative path
    if (cfg.run_mode)
      cfg.tmp_dir = cfg.output_root + cfg.tmp_dir;
    else
      cfg.tmp_dir = cfg.output_dir + cfg.tmp_dir;
  }

  load(cfg.main_iteration, pt, "main_iteration");

  load(cfg.entry_point, pt, "entry_point");

  load(cfg.use_additional_contigs, pt, "use_additional_contigs");
  load(cfg.use_unipaths, pt, "use_unipaths");

  load(cfg.pb, pt, "pacbio_processor");

  load(cfg.additional_contigs, pt, "additional_contigs");

  load(cfg.rr_enable, pt, "rr_enable");
  load(cfg.single_reads_rr, pt, "single_reads_rr");
  cfg.use_single_reads = false;

  load(cfg.mismatch_careful, pt, "mismatch_careful");
  load(cfg.correct_mismatches, pt, "correct_mismatches");
  load(cfg.paired_info_statistics, pt, "paired_info_statistics");
  load(cfg.paired_info_scaffolder, pt, "paired_info_scaffolder");
  load(cfg.cut_bad_connections, pt, "cut_bad_connections");
  load(cfg.gap_closer_enable, pt, "gap_closer_enable");

  load(cfg.max_repeat_length, pt, (cfg.ds.single_cell ? "max_repeat_length_sc" : "max_repeat_length"));

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

  load(cfg.max_threads, pt, "max_threads");
  // Fix number of threads according to OMP capabilities.
  cfg.max_threads = std::min(cfg.max_threads, (size_t) omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads((int) cfg.max_threads);

  load(cfg.max_memory, pt, "max_memory");

  load(cfg.diploid_mode, pt,	"diploid_mode");

  path::CheckFileExistenceFATAL(cfg.dataset_file);
  boost::property_tree::ptree ds_pt;
  boost::property_tree::read_info(cfg.dataset_file, ds_pt);
  load(cfg.ds, ds_pt, true);

  load(cfg.ade, pt, (cfg.ds.single_cell ? "sc_ade" : "usual_ade")); // advanced distance estimator:

  load(cfg.pos, pt, "pos"); // position handler:

  load(cfg.est_mode, pt, "estimation_mode");

  load(cfg.amb_de, pt, "amb_de");
  cfg.amb_de.enabled = (cfg.diploid_mode) ? true : false;

  load(cfg.rm, pt, "resolving_mode");

  if (cfg.rm == rm_path_extend) {
      load(cfg.de, pt, (cfg.ds.single_cell ? "sc_de" : "usual_de"));
  }
  else {
      load(cfg.de, pt, (cfg.ds.single_cell ? "old_sc_de" : "old_usual_de"));
  }

  cfg.pe_params.name = cfg.ds.single_cell ? "singlecell" : "multicell";
  load(cfg.pe_params, pt, "path_extend_params");
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

  load(cfg.con, pt, "construction");
  load(cfg.gc, pt, "gap_closer");
  load(cfg.graph_read_corr, pt, "graph_read_corr");
  load(cfg.kcm, pt, "kmer_coverage_model");
  load(cfg.need_consensus, pt, "need_consensus");
  load(cfg.uncorrected_reads, pt, "uncorrected_reads");
  load(cfg.mismatch_ratio, pt, "mismatch_ratio");

  load(cfg.con, pt, "construction");
  load(cfg.sensitive_map, pt, "sensitive_mapper");
  load(cfg.flanking_range, pt, "flanking_range");

  load(cfg.info_printers, pt, "info_printers");
  if (!cfg.developer_mode) {
      for (auto iter = cfg.info_printers.begin(); iter != cfg.info_printers.end(); ++iter) {
          clear(iter->second);
      }
  }
  load_reads(cfg.ds, cfg.input_dir);
  load_reference_genome(cfg.ds, cfg.input_dir);

  cfg.need_mapping = cfg.developer_mode || cfg.correct_mismatches 
                        || cfg.gap_closer_enable || cfg.rr_enable;

  load(cfg.simp, pt, "default");

  if (cfg.ds.single_cell)
    load(cfg.simp, pt, "sc", false);

  if (cfg.mismatch_careful)
    load(cfg.simp, pt, "careful", false);

  if (cfg.diploid_mode) {
    load(cfg.simp, pt, "diploid_simp", false);
  }

}

void load(debruijn_config& cfg, const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(filename, pt);

  load(cfg, pt, true);
}

};

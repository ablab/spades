#include "pe_config_struct.hpp"
#include "config_common.hpp"

namespace path_extend {

void load(output_broken_scaffolds& obs, boost::property_tree::ptree const& pt, std::string const& key, bool complete) {

  if (complete || pt.find(key) != pt.not_found()) {
    std::string ep = pt.get<std::string>(key);
    obs = pe_config::output_broken_scaffolds_id(ep);
  }

}

void load(pe_config::OutputParamsT& o, boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(o.write_overlaped_paths,   pt, "write_overlaped_paths" );
  load(o.write_paths,             pt, "write_paths"           );
  load(o.write_path_loc,          pt, "write_path_loc"        );
}

void load(pe_config::VisualizeParamsT& o, boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;
  load(o.print_overlaped_paths,   pt, "print_overlaped_paths" );
  load(o.print_paths,             pt, "print_paths"           );
}

void load(pe_config::UtilsT& u, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;
  load(u.mode, pt, "mode");
  load(u.file1, pt, "file1");
  load(u.file2, pt, "file2");

  load(u.advanced, pt, "advanced");
  load(u.clustered, pt, "clustered");
  load(u.insert_size, pt, "insert_size");
  load(u.read_size, pt, "read_size");
  load(u.dev, pt, "dev");
}


void load(pe_config::DatasetT::PairedLibT& pl, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;
  load(pl.read_size  , pt, "read_size"  );
  load(pl.insert_size, pt, "insert_size");
  load(pl.var        , pt, "var"        );
  load(pl.path       , pt, "path"       );
}

void load(pe_config::DatasetT& ds, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;

  //ds.reference_genom = pt.get_optional<std::string>("reference_genome");

  load(ds.graph_file, pt, "graph_file");

  load(ds.param_set,  pt, "param_set");
  //TODO
  //load(ds.libs,       pt,     "libs" );
}


void load(pe_config::ParamSetT::ExtensionOptionsT::SelectOptionsT& so, boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(so.priority_coeff,     pt, "priority_coeff");
  load(so.weight_threshold,   pt, "weight_threshold");
  so.single_threshold = pt.get_optional<double>("single_threshold");
}

void load(pe_config::ParamSetT::ExtensionOptionsT& es, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;

  load(es.try_deep_search , pt, "try_deep_search");
  load(es.select_options  , pt, es.metric.c_str()          );
}

void load(pe_config::ParamSetT::LoopRemovalT& lr, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;
  load(lr.inspect_short_loops, pt,"inspect_short_loops"      );

  load(lr.max_loops          , pt,"max_loops"          );
  load(lr.full_loop_removal  , pt,"full_loop_removal"  );
}


void load(pe_config::ParamSetT::FilterOptionsT& fo, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;
  load(fo.remove_overlaps      , pt, "remove_overlaps"      );
}


void load(pe_config::ParamSetT::ScaffolderOptionsT& so, boost::property_tree::ptree const& pt, bool /*complete*/)
{
  using config_common::load;
  load(so.on      , pt, "on"      );
  load(so.cutoff      , pt, "cutoff"      );
  load(so.rel_cutoff      , pt, "rel_cutoff"      );
  load(so.sum_threshold      , pt, "sum_threshold"      );

  load(so.cluster_info      , pt, "cluster_info"      );
  load(so.cl_threshold      , pt, "cl_threshold"      );

  load(so.fix_gaps      , pt, "fix_gaps"      );
  load(so.min_gap_score      , pt, "min_gap_score"      );
  load(so.max_must_overlap      , pt, "max_must_overlap"      );
  load(so.max_can_overlap      , pt, "max_can_overlap"      );
  load(so.short_overlap      , pt, "short_overlap"      );
  load(so.artificial_gap      , pt, "artificial_gap"      );
}

void load(pe_config::ParamSetT& p, boost::property_tree::ptree const& pt, bool /*complete*/) {

  using config_common::load;

  load(p.metric, pt,  "metric");
  load(p.normalize_weight, pt,  "normalize_weight");
  load(p.normalize_by_coverage, pt,  "normalize_by_coverage");

  load(p.improve_paired_info, pt,  "improve_paired_info");

  load(p.split_edge_length, pt, "split_edge_length");

  p.extension_options.metric = p.metric;
  p.mate_pair_options.metric = "path_cover";

  load(p.extension_options, pt, "extension_options");
  load(p.mate_pair_options, pt, "mate_pair_options");
  load(p.scaffolder_options, pt, "scaffolder");
  load(p.loop_removal,      pt, "loop_removal");
  load(p.filter_options,    pt, "filter_options");
}

void load(pe_config::LongReads& p, boost::property_tree::ptree const& pt, bool /*complete*/) {

  using config_common::load;
  load(p.filtering, pt, "filtering");
  load(p.priority, pt, "priority");
}

void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(p.debug_output, pt,  "debug_output"   );

  load(p.output      , pt,  "output"   );
  load(p.viz         , pt,  "visualize");
  load(p.param_set   , pt,  p.name.c_str()   );
  load(p.obs         , pt,  "output_broken_scaffolds");
  load(p.long_reads, pt, "long_reads");


  if (!p.debug_output) {
    p.output.DisableAll();
    p.viz.DisableAll();
  }
  p.etc_dir = "path_extend";
}


// main long contigs config load function
void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool /*complete*/) {
  using config_common::load;

  load(pe_cfg.dataset_name           , pt, "dataset"               );
  load(pe_cfg.dataset                , pt,  pe_cfg.dataset_name.c_str()  );

  pe_cfg.params.name = pe_cfg.dataset.param_set;
  load(pe_cfg.params                 , pt, "pe_params"             );
}

};


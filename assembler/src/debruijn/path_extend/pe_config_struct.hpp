//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * lc_config_struct.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef LC_CONFIG_STRUCT_HPP_
#define LC_CONFIG_STRUCT_HPP_

#include "config_singl.hpp"
#include "cpp_utils.hpp"

#include <boost/optional.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/bimap.hpp>

#include <string>
#include <vector>

namespace path_extend {

enum output_broken_scaffolds {
    obs_none,
    obs_break_gaps,
    obs_break_all
};

enum scaffolding_mode {
    sm_old,
    sm_2015,
    sm_combined,
    sm_old_pe_2015
};

inline bool is_2015_scaffolder_enabled(const scaffolding_mode mode) {
    return (mode != sm_old);
}

// struct for path extend subproject's configuration file
struct pe_config {

  typedef boost::bimap<std::string, output_broken_scaffolds> output_broken_scaffolds_id_mapping;

  static const output_broken_scaffolds_id_mapping FillOBSInfo() {
    output_broken_scaffolds_id_mapping::value_type info[] = {
              output_broken_scaffolds_id_mapping::value_type("none", obs_none),
              output_broken_scaffolds_id_mapping::value_type("break_gaps", obs_break_gaps),
              output_broken_scaffolds_id_mapping::value_type("break_all", obs_break_all)
    };

    return output_broken_scaffolds_id_mapping(info, utils::array_end(info));
  }

  static const output_broken_scaffolds_id_mapping& output_broken_scaffolds_info() {
    static output_broken_scaffolds_id_mapping output_broken_scaffolds_info = FillOBSInfo();
    return output_broken_scaffolds_info;
  }

  static const std::string& output_broken_scaffolds_name(output_broken_scaffolds obs) {
    auto it = output_broken_scaffolds_info().right.find(obs);
    VERIFY_MSG(it != output_broken_scaffolds_info().right.end(),
               "No name for output broken scaffolds mode id = " << obs);

    return it->second;
  }

  static output_broken_scaffolds output_broken_scaffolds_id(std::string name) {
    auto it = output_broken_scaffolds_info().left.find(name);
    VERIFY_MSG(it != output_broken_scaffolds_info().left.end(),
               "There is no output broken scaffolds mode with name = " << name);

    return it->second;
  }

  typedef boost::bimap<std::string, scaffolding_mode> scaffolding_mode_id_mapping;

  static const scaffolding_mode_id_mapping FillSMInfo() {
      scaffolding_mode_id_mapping::value_type info[] = {
              scaffolding_mode_id_mapping::value_type("old", sm_old),
              scaffolding_mode_id_mapping::value_type("2015", sm_2015),
              scaffolding_mode_id_mapping::value_type("combined", sm_combined),
              scaffolding_mode_id_mapping::value_type("old_pe_2015", sm_old_pe_2015)
    };

    return scaffolding_mode_id_mapping(info, utils::array_end(info));
  }

  static const scaffolding_mode_id_mapping& scaffolding_mode_info() {
    static scaffolding_mode_id_mapping scaffolding_mode_info = FillSMInfo();
    return scaffolding_mode_info;
  }

  static const std::string& scaffolding_mode_name(scaffolding_mode sm) {
    auto it = scaffolding_mode_info().right.find(sm);
    VERIFY_MSG(it != scaffolding_mode_info().right.end(),
               "No name for scaffolding mode id = " << sm);

    return it->second;
  }

  static scaffolding_mode scaffolding_mode_id(std::string name) {
    auto it = scaffolding_mode_info().left.find(name);
    VERIFY_MSG(it != scaffolding_mode_info().left.end(),
               "There is no scaffolding mode with name = " << name);

    return it->second;
  }

  struct OutputParamsT {
    bool write_overlaped_paths;
    bool write_paths;

    void DisableAll() {
      write_overlaped_paths = false;
      write_paths = false;
    }
  };



  struct VisualizeParamsT {
    bool print_overlaped_paths;
    bool print_paths;

    void DisableAll() {
      print_overlaped_paths = false;
      print_paths = false;
    }
  };

  struct ParamSetT {
    scaffolding_mode sm;

    bool normalize_weight;
    size_t split_edge_length;
    bool cut_all_overlaps;

    struct ExtensionOptionsT {
        bool use_default_single_threshold;
        double single_threshold;
        double weight_threshold;
        double priority_coeff;
    } extension_options;

    ExtensionOptionsT mate_pair_options;


    struct ScaffolderOptionsT {
      bool on;
      int cutoff;
      double rel_cutoff;
      double sum_threshold;

      bool cluster_info;
      double cl_threshold;

      bool fix_gaps;
      double min_gap_score;
      double max_must_overlap;
      double max_can_overlap;
      int short_overlap;
      size_t artificial_gap;

      bool use_old_score;

      size_t min_overlap_length;
      double flank_addition_coefficient;
      double flank_multiplication_coefficient;
    } scaffolder_options;


    struct LoopRemovalT {
      size_t max_loops;
      size_t mp_max_loops;
    } loop_removal;

    bool remove_overlaps;
    bool use_coordinated_coverage;

    struct CoordinatedCoverageT {
      size_t max_edge_length_in_repeat;
      double delta;
    } coordinated_coverage;
      struct Scaffolding2015 {
          bool autodetect;
          size_t min_unique_length;
          double unique_coverage_variation;
      } scaffolding2015;
      struct ScaffoldGraphParamsT {
          bool construct;
          bool output;
          size_t min_read_count;
          bool graph_connectivity;
          size_t max_path_length;
      } scaffold_graph_params;
  };

  struct LongReads {
	  double filtering;
	  double weight_priority;
	  double unique_edge_priority;
  };

  struct AllLongReads{
      LongReads single_reads;
      LongReads pacbio_reads;
      LongReads contigs;
  };


  struct MainPEParamsT {
    output_broken_scaffolds obs;

    bool finalize_paths;
    bool debug_output;
    std::string etc_dir;

    OutputParamsT output;
    VisualizeParamsT viz;
    ParamSetT param_set;
    AllLongReads long_reads;
  }; // params;


//  std::string dataset_name;


};

void load(pe_config::ParamSetT& p, boost::property_tree::ptree const& pt, bool complete = true);
void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt, bool complete = true);
//void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool complete);

}

//typedef config_common::config<path_extend::pe_config> pe_cfg;

#endif /* CONFIG_STRUCT_HPP_ */

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

const char * const pe_cfg_filename = "./config/debruijn/path_extend/lc_config.info";

enum output_broken_scaffolds {
    obs_none,
    obs_break_gaps,
    obs_break_all
};

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
               "No name for working stage id = " << obs);

    return it->second;
  }

  static output_broken_scaffolds output_broken_scaffolds_id(std::string name) {
    auto it = output_broken_scaffolds_info().left.find(name);
    VERIFY_MSG(it != output_broken_scaffolds_info().left.end(),
               "There is no working stage with name = " << name);

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
    bool normalize_weight;
    size_t split_edge_length;

    struct ExtensionOptionsT {
        bool recalculate_threshold;
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
      int artificial_gap;
    } scaffolder_options;


    struct LoopRemovalT {
      size_t max_loops;
      size_t mp_max_loops;
    } loop_removal;

    bool remove_overlaps;
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
    std::string name;
    output_broken_scaffolds obs;

    bool debug_output;
    std::string etc_dir;

    OutputParamsT output;
    VisualizeParamsT viz;
    ParamSetT param_set;
    AllLongReads long_reads;
    bool cut_all_overlaps;
  } params;

  std::string dataset_name;


};



void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt, bool complete);
void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool complete);

}

typedef config_common::config<path_extend::pe_config> pe_cfg;

#endif /* CONFIG_STRUCT_HPP_ */

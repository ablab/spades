//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

#include <boost/optional.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include <string>
#include <vector>

namespace path_extend {

const char * const pe_cfg_filename = "./config/debruijn/path_extend/lc_config.info";

// struct for long_contigs subproject's configuration file
struct pe_config {
  struct DatasetT {
    struct PairedLibT {
      size_t read_size;
      size_t insert_size;
      size_t var;

      std::string path;
    };

    std::string param_set;
    std::string graph_file;

    std::vector<PairedLibT> libs;
    //boost::optional<std::string> reference_genome;
  };

  struct OutputParamsT {
    bool write_seeds;
    bool write_overlaped_paths;
    bool write_paths;
    bool write_path_loc;

    void DisableAll() {
      write_seeds = false;
      write_overlaped_paths = false;
      write_paths = false;
      write_path_loc = false;
    }
  };

  struct VisualizeParamsT {
    bool print_seeds;
    bool print_overlaped_paths;
    bool print_paths;

    void DisableAll() {
      print_seeds = false;
      print_overlaped_paths = false;
      print_paths = false;
    }
  };

  struct ResearchT {
    bool on;

    bool count_seed_weight;
    bool count_path_weight;

    bool fiter_seeds_by_id;
    std::vector<size_t> seed_ids;
  };

  struct ParamSetT {
    std::string metric;
    bool normalize_weight;
    bool normalize_by_coverage;

    bool improve_paired_info;

    struct SeedSelectionT {
      std::string metric;

      double min_coverage;
      double start_egde_coverage;
      size_t max_cycles;

      bool exclude_chimeric;
      int chimeric_delta;

      bool check_trusted;
      double threshold;
    } seed_selection;


    struct ExtensionOptionsT {
      std::string metric;

      bool try_deep_search;

      struct SelectOptionsT {
        boost::optional<double> single_threshold;
        double weight_threshold;
        double priority_coeff;

        SelectOptionsT() {}
        SelectOptionsT(const SelectOptionsT& so)
            : single_threshold(so.single_threshold), weight_threshold(so.weight_threshold), priority_coeff(so.priority_coeff) {}
      } select_options;

    } extension_options;


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
      bool inspect_short_loops;

      size_t max_loops;
      bool full_loop_removal;
    } loop_removal;


    struct FilterOptionsT {
      bool remove_overlaps;
    } filter_options;
  };

  struct UtilsT {
    int mode;
    std::string file1;
    std::string file2;

    std::string clustered;
    std::string advanced;
    size_t insert_size;
    size_t read_size;
    size_t dev;
  };

  struct MainPEParamsT {
    std::string name;

    bool debug_output;

    OutputParamsT output;
    VisualizeParamsT viz;
    ParamSetT param_set;
  } params;

  std::string dataset_name;
  DatasetT dataset;
};

void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt, bool complete);
void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool complete);
}

typedef config_common::config<path_extend::pe_config> pe_cfg;

#endif /* CONFIG_STRUCT_HPP_ */

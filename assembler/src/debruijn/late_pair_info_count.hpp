//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "simplification.hpp"
#include "graph_construction.hpp"
#include "dataset_readers.hpp"

#include "de/insert_size_refiner.hpp"
#include "de/paired_info.hpp"

namespace debruijn_graph {

  void check_new_index(conj_graph_pack& gp);

  void late_pair_info_count(conj_graph_pack& gp, PairedIndexT& paired_index)
  {
    exec_simplification(gp);

    if (!cfg::get().developer_mode) {
      paired_index.Attach();
      paired_index.Init();
    }

    EdgeQuality<Graph> edge_qual(gp.g, gp.index,
        gp.kmer_mapper, gp.genome);
    total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
    total_labeler tot_lab(&graph_struct);
    CompositeLabeler<Graph> labeler(tot_lab, edge_qual);

    auto mapper = MapperInstance(gp);

    ChimericEdgesLengthStats<Graph> chimeric_stats(gp.g, edge_qual, mapper->MapSequence(gp.genome), 100000);

    chimeric_stats();

    if (cfg::get().paired_mode) {
      typedef io::ReadStreamVector<SequencePairedReadStream> MultiStreamType;
      typedef io::ReadStreamVector<PairedReadStream> SingleStreamType;
      size_t edge_length_threshold = Nx(gp.g, 50);

      INFO("STAGE == Counting Late Pair Info");

      if (cfg::get().use_multithreading) {
        auto streams = paired_binary_readers(false, 0);
        bool success = RefineInsertSize(gp, *streams, cfg::get_writable(), edge_length_threshold);
        if (!success)
          return;

        auto paired_streams = paired_binary_readers(true, *cfg::get().ds.IS);
        FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
            gp.kmer_mapper, paired_index, *paired_streams, gp.k_value);
      } else {
        auto_ptr<PairedReadStream> stream = paired_easy_reader(false, 0);
        SingleStreamType streams(stream.get());
        bool success = RefineInsertSize(gp, streams, cfg::get_writable(), edge_length_threshold);
        if (!success)
          return;

        auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(true, *cfg::get().ds.IS);
        SingleStreamType paired_streams(paired_stream.get());

        FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
            gp.kmer_mapper, paired_index, paired_streams, gp.k_value);
      }
    }
  }

  void load_late_pair_info_count(conj_graph_pack& gp,
      PairedIndexT& paired_index, path::files_t* used_files)
  {
    string p = path::append_path(cfg::get().load_from, "late_pair_info_counted");
    used_files->push_back(p);

    ScanWithPairedIndex(p, gp, paired_index);
    load_estimated_params(p);
  }

  void save_late_pair_info_count(conj_graph_pack& gp, PairedIndexT& paired_index) {
    if (cfg::get().make_saves || (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles && cfg::get().paired_mode)) {
      if (!cfg::get().make_saves)
        make_dir(cfg::get().output_saves);

      string p = path::append_path(cfg::get().output_saves, "late_pair_info_counted");
      INFO("Saving current state to " << p);
      PrintWithPairedIndex(p, gp, paired_index);
      write_estimated_params(p);
    }
  }

  void exec_late_pair_info_count(conj_graph_pack& gp, PairedIndexT& paired_index) {
    if (cfg::get().entry_point <= ws_late_pair_info_count) {
      late_pair_info_count(gp, paired_index);
      save_late_pair_info_count(gp, paired_index);
    } else {
      INFO("Loading Late Pair Info Count");
      path::files_t used_files;
      load_late_pair_info_count(gp, paired_index, &used_files);
      link_files_by_prefix(used_files, cfg::get().output_saves);
    }
  }

}

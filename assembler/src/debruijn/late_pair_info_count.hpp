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

  typedef io::ReadStreamVector<SequencePairedReadStream> MultiStreamType;
  typedef io::ReadStreamVector<PairedReadStream> SingleStreamType;

  void check_new_index(conj_graph_pack& gp);

  void late_pair_info_count(conj_graph_pack& gp, PairedIndexT& paired_index)
  {
    exec_simplification(gp);

	/*for ( auto index_iterator = gp.index.inner_index().value_begin(); index_iterator < gp.index.inner_index().value_end(); ++index_iterator ){
              std::cout << index_iterator->edgeId_ << " " << index_iterator->offset_ << " " << index_iterator->count_ << std::endl;
     }*/


    if (!cfg::get().developer_mode) {
      paired_index.Attach();
      paired_index.Init();
    }

    if (cfg::get().paired_mode) {

      size_t edge_length_threshold = Nx(gp.g, 50);

      INFO("STAGE == Counting Late Pair Info");

      if (cfg::get().use_multithreading) {
        auto streams = paired_binary_readers(false, 0);
        bool success = RefineInsertSize(gp, *streams, cfg::get_writable(), edge_length_threshold);
        if (!success)
          return;

        auto paired_streams = paired_binary_readers(true, cfg::get().ds.IS());
        FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
                                           gp.kmer_mapper, paired_index, *paired_streams, gp.k_value);
      } else {
        auto_ptr<PairedReadStream> stream = paired_easy_reader(false, 0);
        SingleStreamType streams(stream.get());
        stream.release();
        bool success = RefineInsertSize(gp, streams, cfg::get_writable(), edge_length_threshold);
        if (!success)
          return;

        auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(true, cfg::get().ds.IS());
        SingleStreamType paired_streams(paired_stream.get());
        paired_stream.release();

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


  // ===== Multiple lib version =====

  void late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices)
  {
    exec_simplification(gp);

    if (!cfg::get().developer_mode) {
        paired_indices.Attach();
        paired_indices.Init();
    }

    if (cfg::get().paired_mode) {
        size_t edge_length_threshold = Nx(gp.g, 50);
        INFO("STAGE == Counting Late Pair Info");

        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
           if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||
                   cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {

               bool success;
               if (cfg::get().use_multithreading) {
                   auto streams = paired_binary_readers(cfg::get().ds.reads[i], false, 0);
                   success = RefineInsertSizeForLib(gp, *streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
               }
               else {
                   auto_ptr<PairedReadStream> stream = paired_easy_reader(cfg::get().ds.reads[i], false, 0);
                   SingleStreamType streams(stream.get());
                   streams.release();
                   success = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
               }

               if (!success) {
                   INFO("Unable to estimate insert size for paired library #" << i);
                   continue;
               }
               else {
                   INFO("Estimated insert size for paired library #" << i);
                   INFO("Insert size = " << cfg::get().ds.reads[i].data().mean_insert_size << ", deviation = " << cfg::get().ds.reads[i].data().insert_size_deviation)
               }

               if (cfg::get().use_multithreading) {
                   auto paired_streams = paired_binary_readers(cfg::get().ds.reads[i], true, cfg::get().ds.reads[i].data().mean_insert_size);
                   FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
                                                  gp.kmer_mapper, paired_indices[i], *paired_streams, gp.k_value);
               }
               else {
                   auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(cfg::get().ds.reads[i], true, cfg::get().ds.reads[i].data().mean_insert_size);
                   SingleStreamType paired_streams(paired_stream.get());
                   paired_stream.release();
                   FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
                                                      gp.kmer_mapper, paired_indices[i], paired_streams, gp.k_value);
               }

           }
        }

    }
  }


  void load_late_pair_info_count(conj_graph_pack& gp,
          PairedIndicesT& paired_indices, path::files_t* used_files)
  {
    string p = path::append_path(cfg::get().load_from, "late_pair_info_counted");
    used_files->push_back(p);

    ScanWithPairedIndices(p, gp, paired_indices);
    load_lib_data(p);
  }

  void save_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices) {
    if (cfg::get().make_saves || (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles && cfg::get().paired_mode)) {
      if (!cfg::get().make_saves)
        make_dir(cfg::get().output_saves);

      string p = path::append_path(cfg::get().output_saves, "late_pair_info_counted");
      INFO("Saving current state to " << p);

      PrintWithPairedIndices(p, gp, paired_indices);
      write_lib_data(p);
    }
    // for informing spades.py about estimated params
    write_lib_data(cfg::get().output_dir + "/");
  }

  void exec_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices) {
    if (cfg::get().entry_point <= ws_late_pair_info_count) {
      late_pair_info_count(gp, paired_indices);
      save_late_pair_info_count(gp, paired_indices);
    } else {
      INFO("Loading Late Pair Info Count");
      path::files_t used_files;
      load_late_pair_info_count(gp, paired_indices, &used_files);
      link_files_by_prefix(used_files, cfg::get().output_saves);
    }
  }


}

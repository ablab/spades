//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * construction.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "io/easy_reader.hpp"
#include "io/vector_reader.hpp"
#include "omni_labelers.hpp"
#include "dataset_readers.hpp"
#include "positions.hpp"

namespace debruijn_graph {

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
    total_labeler& tl, PairedIndexT& paired_index);

} // namespace debruijn_graph

// todo: move impl to *.cpp

namespace debruijn_graph {

template<class Read>
void construct_graph(io::ReadStreamVector< io::IReader<Read> >& streams,
                     conj_graph_pack& gp, ReadStream* contigs_stream = 0) {
  INFO("STAGE == Constructing Graph");

  debruijn_config::construction params = cfg::get().con;
  params.early_tc.enable &= !cfg::get().ds.single_cell;

  //  size_t rl = ConstructGraphWithCoverage(cfg::get().K, streams, gp.g, gp.index, contigs_stream);
  size_t rl = ConstructGraphWithCoverage(cfg::get().K, params, streams, gp.g,
      gp.index, contigs_stream);

  if (!cfg::get().ds.RL()) {
    INFO("Figured out: read length = " << rl);
    cfg::get_writable().ds.set_RL(rl);
  } else if (cfg::get().ds.RL() != rl)
    WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL() << ", not " << rl);
}

void load_construction(conj_graph_pack& gp, path::files_t* files) {
  string p = path::append_path(cfg::get().load_from, "constructed_graph");
  files->push_back(p);
  ScanGraphPack(p, gp);
  load_lib_data(p);
}

void save_construction(conj_graph_pack& gp) {
  if (cfg::get().make_saves) {
    string p = path::append_path(cfg::get().output_saves, "constructed_graph");
    INFO("Saving current state to " << p);
    PrintGraphPack(p, gp);
    write_lib_data(p);
  }
}

void exec_construction(conj_graph_pack& gp) {
  if (cfg::get().entry_point <= ws_construction) {

//    if (cfg::get().etalon_graph_mode) {
//      typedef io::VectorReader<io::SingleRead> GenomeStream;
//      GenomeStream genome_stream(io::SingleRead("genome", gp.genome.str()));
//      std::vector <ReadStream*> streams(1, &genome_stream);
//      construct_graph(streams, gp);
//    } else

    //has to be separate stream for not counting it in coverage
    ReadStream* additional_contigs_stream = 0;
    if (cfg::get().use_additional_contigs) {
      INFO("Contigs from previous K will be used");
      additional_contigs_stream = new io::EasyReader(
          cfg::get().additional_contigs, true);
    }

    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd || cfg::get().ds.reads[i].type() == io::LibraryType::SingleReads) {
            libs_for_construction.push_back(i);
        }
    }

    if (cfg::get().use_multithreading) {
      auto streams = single_binary_readers_for_libs(libs_for_construction, true, true);
      construct_graph<io::SingleReadSeq>(*streams, gp, additional_contigs_stream);
    } else {
      auto single_stream = single_easy_reader_for_libs(libs_for_construction, true, true);
      io::ReadStreamVector<ReadStream> streams(single_stream.get());
      single_stream.release();
      construct_graph<io::SingleRead>(streams, gp, additional_contigs_stream);
    }

    save_construction(gp);
  } else {
    INFO("Loading Construction");

        path::files_t used_files;
    load_construction(gp, &used_files);
    link_files_by_prefix(used_files, cfg::get().output_saves);
//    OnlineVisualizer online(gp);
//    online.run();
  }

  if (cfg::get().developer_mode) {
    if (gp.genome.size() > 0) {
      FillPos(gp, gp.genome, "ref0");
      FillPos(gp, !gp.genome, "ref1");
    }

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

} //namespace debruijn_graph

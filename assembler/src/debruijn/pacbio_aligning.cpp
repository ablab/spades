//****************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "graph_simplification.hpp"
#include "omni/visualization/graph_labeler.hpp"
#include "io/single_read.hpp"
#include <algorithm>
#include "pacbio/pac_index.hpp"
#include "pacbio/pacbio_gap_closer.hpp"
#include "long_read_storage.hpp"

#include "pacbio_aligning.hpp"

namespace debruijn_graph {

void ProcessReadsBatch(conj_graph_pack &gp,
                       std::vector<io::SingleRead>& reads,
                       pacbio::PacBioMappingIndex<ConjugateDeBruijnGraph>& pac_index,
                       PathStorage<Graph>& long_reads, pacbio::GapStorage<Graph>& gaps,
                       size_t buf_size, int n, size_t min_gap_quantity, pacbio::StatsCounter stats) {
    vector<PathStorage<Graph> > long_reads_by_thread(cfg::get().max_threads,
                                                     PathStorage<Graph>(gp.g));
    vector<pacbio::GapStorage<Graph> > gaps_by_thread(cfg::get().max_threads,
                                              pacbio::GapStorage<Graph>(gp.g, min_gap_quantity));
    vector<pacbio::StatsCounter> stats_by_thread(cfg::get().max_threads);


#   pragma omp parallel for shared(reads, long_reads_by_thread, pac_index, n)
    for (size_t i = 0; i < buf_size; ++i) {
        if (i % 1000 == 0) {
            DEBUG("thread number " << omp_get_thread_num());
        }
        size_t thread_num = omp_get_thread_num();
        Sequence seq(reads[i].sequence());
#       pragma omp atomic
        n++;
        auto current_read_mapping = pac_index.GetReadAlignment(seq);
        auto aligned_edges = current_read_mapping.main_storage;
        auto gaps = current_read_mapping.gaps;
        for (auto iter = gaps.begin(); iter != gaps.end(); ++iter)
            gaps_by_thread[thread_num].AddGap(*iter, true);

        for (auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter)
            long_reads_by_thread[thread_num].AddPath(*iter, 1, true);
        //counting stats:
        for (auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter)
            stats_by_thread[thread_num].path_len_in_edges[iter->size()] ++;

#       pragma omp critical
        {
            VERBOSE_POWER(n, " reads processed");
        }
    }

    for (size_t i = 0; i < cfg::get().max_threads; i++) {
        long_reads.AddStorage(long_reads_by_thread[i]);
        gaps.AddStorage(gaps_by_thread[i]);
        stats.AddStorage(stats_by_thread[i]);
    }
}

void align_pacbio(conj_graph_pack &gp, int lib_id) {
    auto pacbio_read_stream = single_easy_reader(cfg::get().ds.reads[lib_id],
                                                 false, false);
    io::ReadStreamList<io::SingleRead> streams(pacbio_read_stream);
 //   pacbio_read_stream.release();
    int n = 0;
    PathStorage<Graph>& long_reads = gp.single_long_reads[lib_id];
    pacbio::StatsCounter stats;
    size_t min_gap_quantity = 2;
    if (cfg::get().ds.reads[lib_id].type() == io::LibraryType::PacBioReads || cfg::get().ds.reads[lib_id].type() == io::LibraryType::SangerReads) {
        min_gap_quantity = cfg::get().pb.pacbio_min_gap_quantity;
    } else {
        min_gap_quantity = cfg::get().pb.contigs_min_gap_quantity;
    }
    pacbio::GapStorage<ConjugateDeBruijnGraph> gaps(gp.g, min_gap_quantity);
    size_t read_buffer_size = 50000;
    std::vector<io::SingleRead> reads(read_buffer_size);
    io::SingleRead read;
    size_t buffer_no = 0;
    INFO("Usign seed size: " << cfg::get().pb.pacbio_k);
    pacbio::PacBioMappingIndex<ConjugateDeBruijnGraph> pac_index(gp.g,
                                                         cfg::get().pb.pacbio_k,
                                                         cfg::get().K, cfg::get().pb.ignore_middle_alignment);
//    path_extend::ContigWriter cw(gp.g);
//    cw.writeEdges("before_rr_with_ids.fasta");
//    ofstream filestr("pacbio_mapped.mpr");
//    filestr.close();
    for (auto iter = streams.begin(); iter != streams.end(); ++iter) {
        auto &stream = *iter;
        while (!stream.eof()) {
            size_t buf_size = 0;
            for (; buf_size < read_buffer_size && !stream.eof(); ++buf_size)
                stream >> reads[buf_size];
            INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");
            DEBUG("master thread number " << omp_get_thread_num());
            ProcessReadsBatch(gp, reads, pac_index, long_reads, gaps, buf_size, n, min_gap_quantity, stats);
            INFO("Processed batch " << buffer_no);
            ++buffer_no;
        }
    }
    map<EdgeId, EdgeId> replacement;
    long_reads.DumpToFile(cfg::get().output_saves + "long_reads_before_rep.mpr",
                          replacement);
    gaps.DumpToFile(cfg::get().output_saves + "gaps.mpr");
    gaps.PadGapStrings();
    gaps.DumpToFile(cfg::get().output_saves +  "gaps_padded.mpr");
    pacbio::PacbioGapCloser<Graph> gap_closer(gp.g);
    gap_closer.ConstructConsensus(cfg::get().max_threads, gaps);
    gap_closer.CloseGapsInGraph(replacement);
    long_reads.ReplaceEdges(replacement);
    for(int j = 0; j < lib_id; j++) {
        gp.single_long_reads[j].ReplaceEdges(replacement);
    }

    gap_closer.DumpToFile(cfg::get().output_saves + "gaps_pb_closed.fasta");
    INFO("Index refill");
    gp.index.Refill();
    INFO("Index refill after PacBio finished");
    if (!gp.index.IsAttached())
        gp.index.Attach();
    INFO("PacBio test finished");
    return;
}

void PacBioAligning::run(conj_graph_pack &gp, const char*) {
    using namespace omnigraph;
    omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    int lib_id = -1;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        io::LibraryType type = cfg::get().ds.reads[i].type();
        if (type == io::LibraryType::PacBioReads ||
                type == io::LibraryType::SangerReads ||
                type == io::LibraryType::TrustedContigs ||
                type == io::LibraryType::UntrustedContigs) {
            lib_id = (int) i;
            align_pacbio(gp, lib_id);
        }
    }

    if (lib_id == -1)
        INFO("no PacBio lib found");

    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(ipp_final_gap_closed);
}

}


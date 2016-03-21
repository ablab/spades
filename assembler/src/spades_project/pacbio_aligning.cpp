//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "data_structures/assembly_graph/graph_alignment/pacbio/pac_index.hpp"
#include "data_structures/assembly_graph/graph_alignment/pacbio/pacbio_gap_closer.hpp"
#include "data_structures/assembly_graph/graph_alignment/long_read_storage.hpp"
#include "io/reads_io/wrapper_collection.hpp"
#include "data_structures/assembly_graph/graph_support/stats/picture_dump.hpp"
#include "pacbio_aligning.hpp"

namespace debruijn_graph {

void ProcessReadsBatch(conj_graph_pack &gp,
                       std::vector<io::SingleRead>& reads,
                       pacbio::PacBioMappingIndex<ConjugateDeBruijnGraph>& pac_index,
                       PathStorage<Graph>& long_reads, pacbio::GapStorage<Graph>& gaps,
                       size_t buf_size, int n, size_t min_gap_quantity, pacbio::StatsCounter& stats) {
    vector<PathStorage<Graph> > long_reads_by_thread(cfg::get().max_threads,
                                                     PathStorage<Graph>(gp.g));
    vector<pacbio::GapStorage<Graph> > gaps_by_thread(cfg::get().max_threads,
                                              pacbio::GapStorage<Graph>(gp.g, min_gap_quantity));
    vector<pacbio::StatsCounter> stats_by_thread(cfg::get().max_threads);

    size_t longer_500 = 0;
    size_t aligned = 0;
    size_t nontrivial_aligned = 0;

#   pragma omp parallel for shared(reads, long_reads_by_thread, pac_index, n, aligned, nontrivial_aligned)
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
        for (auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter) {
            stats_by_thread[thread_num].path_len_in_edges[iter->size()]++;
        }
#       pragma omp critical
        {
//            INFO(current_read_mapping.seed_num);
            if (seq.size() > 500) {
                longer_500++;
                if (aligned_edges.size() > 0) {
                    aligned++;
                    stats_by_thread[thread_num].seeds_percentage[size_t(
                            floor(double(current_read_mapping.seed_num) * 1000.0 / (double) seq.size()))]++;
                    for (size_t j = 0; j < aligned_edges.size(); j++) {
                        if (aligned_edges[j].size() > 1) {
                            nontrivial_aligned++;
                            break;
                        }
                    }
                }
            }
        }
#       pragma omp critical
        {
            VERBOSE_POWER(n, " reads processed");
        }
    }
    INFO("Read batch of size: " << buf_size << " processed; "<< longer_500 << " of them longer than 500; among long reads aligned: " << aligned << "; paths of more than one edge received: " << nontrivial_aligned );

    for (size_t i = 0; i < cfg::get().max_threads; i++) {
        long_reads.AddStorage(long_reads_by_thread[i]);
        gaps.AddStorage(gaps_by_thread[i]);
        stats.AddStorage(stats_by_thread[i]);
    }
}

void align_pacbio(conj_graph_pack &gp, int lib_id, bool make_additional_saves) {
    io::ReadStreamList<io::SingleRead> streams;
    for (const auto& reads : cfg::get().ds.reads[lib_id].single_reads())
      //do we need input_file function here?
      streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(reads)));

    //make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(file));
    //    auto pacbio_read_stream = single_easy_reader(cfg::get().ds.reads[lib_id],
//    false, false);

//    io::ReadStreamList<io::SingleRead> streams(pacbio_read_stream);
 //   pacbio_read_stream.release();
    int n = 0;
    PathStorage<Graph>& long_reads = gp.single_long_reads[lib_id];
    pacbio::StatsCounter stats;
    size_t min_gap_quantity = 2;
    size_t rtype = 0;
    bool consensus_gap_closing = false;
    if (cfg::get().ds.reads[lib_id].type() == io::LibraryType::PacBioReads || 
        cfg::get().ds.reads[lib_id].type() == io::LibraryType::SangerReads || 
        cfg::get().ds.reads[lib_id].type() == io::LibraryType::NanoporeReads) {
        min_gap_quantity = cfg::get().pb.pacbio_min_gap_quantity;
        rtype = 1;
        consensus_gap_closing = true;
    } else {
        min_gap_quantity = cfg::get().pb.contigs_min_gap_quantity;
        rtype = 2;
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
//    cw.WriteEdges("before_rr_with_ids.fasta");
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
     //       INFO("Processed batch " << buffer_no);
            ++buffer_no;
        }
    }
    string ss = (rtype == 1 ? "long reads": "contigs");
    INFO("For lib " << lib_id << " of " << ss <<" :");
    stats.report();
    map<EdgeId, EdgeId> replacement;
    size_t min_stats_cutoff =(rtype == 1 ? 1  : 0);
    if (make_additional_saves)
        long_reads.DumpToFile(cfg::get().output_saves + "long_reads_before_rep.mpr",
                          replacement, min_stats_cutoff, true);
    gaps.DumpToFile(cfg::get().output_saves + "gaps.mpr");
    gaps.PadGapStrings();
    if (make_additional_saves)
        gaps.DumpToFile(cfg::get().output_saves +  "gaps_padded.mpr");
    pacbio::PacbioGapCloser<Graph> gap_closer(gp.g, consensus_gap_closing);
    gap_closer.ConstructConsensus(cfg::get().max_threads, gaps);
    gap_closer.CloseGapsInGraph(replacement);
    long_reads.ReplaceEdges(replacement);
    for(int j = 0; j < lib_id; j++) {
        gp.single_long_reads[j].ReplaceEdges(replacement);
    }

    gap_closer.DumpToFile(cfg::get().output_saves + "gaps_pb_closed.fasta");
    INFO("PacBio aligning finished");
    return;
}

void PacBioAligning::run(conj_graph_pack &gp, const char*) {
    using namespace omnigraph;
    omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    int lib_id = -1;
    bool make_additional_saves = parent_->saves_policy().make_saves_;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if ( cfg::get().ds.reads[i].is_pacbio_alignable() ) {
            lib_id = (int) i;
            align_pacbio(gp, lib_id, make_additional_saves);
        }
    }

    if (lib_id == -1)
        INFO("no PacBio lib found");

    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(ipp_final_gap_closed);
}

}


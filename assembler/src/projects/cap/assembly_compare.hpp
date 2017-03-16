//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "pipeline/graphio.hpp"
#include "utils/stl_utils.hpp"
#include "modules/simplification/cleaner.hpp"
#include "io/reads/splitting_wrapper.hpp"
#include "io/reads/multifile_reader.hpp"
#include <boost/algorithm/string/predicate.hpp>

#include "coloring.hpp"
#include "colored_graph_construction.hpp"
#include "untangling.hpp"
#include "assembly_problem_detection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "simple_indel_finder.hpp"
#include "test_utils.hpp"

namespace cap {

//class RCSplittingStream: public io::DelegatingReaderWrapper<io::SingleReaGapsRemoverd> {
//private:
//    io::SplittingWrapper filtered_reader_;
//    io::RCReaderWrapper<io::SingleRead> stream_;
//public:
//    RCSplittingStream(ContigStream &base_stream) :
//            filtered_reader_(base_stream), stream_(filtered_reader_) {
//        Init(stream_);
//    }
//};

//todo finish later
//template<class gp_t1, class gp_t2>
//void ConvertToBPGraphPack(const gp_t1& gp
//        , const ColorHandler<typename gp_t1::graph_t>& coloring
//        , gp_t2& bp_gp) {
//    string tmp_dir = "/home/snurk/tmp/";
//    string filename = tmp_dir + "tmp";
//    make_dir(tmp_dir);
//    PrintGraphPack(filename, gp);
//    typename ScannerTraits<typename gp_t2::graph_t>::Scanner scanner(bp_gp.g);
//    ScanBasicGraph(filename, scanner);
//    scanner.loadPositions(filename, bp_gp.edge_pos);
//    //
//}

template<class Graph>
void DeleteEdgesByColor(Graph& g, const ColorHandler<Graph>& coloring,
        TColorSet color) {
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        if (coloring.Color(*it) == color) {
            g.DeleteEdge(*it);
        }
    }
    omnigraph::Cleaner<Graph>(g).Run();
}

template<class Graph>
class GapsRemover {
    typedef typename Graph::VertexId VertexId;
    Graph& g_;
    const ColorHandler<Graph>& coloring_;
    const TColorSet gap_color_;
    const size_t length_bound_;
public:
    GapsRemover(Graph& g, const ColorHandler<Graph>& coloring,
            TColorSet gap_color, size_t length_bound) :
            g_(g), coloring_(coloring), gap_color_(gap_color), length_bound_(
                    length_bound) {

    }

    void RemoveGaps() {
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == gap_color_
                    && g_.length(*it) <= length_bound_
                    && g_.CanCompressVertex(g_.EdgeStart(*it))
                    && g_.CanCompressVertex(g_.EdgeEnd(*it))) {
                VertexId start = g_.EdgeStart(*it);
                VertexId end = g_.EdgeEnd(*it);
                if (!g_.RelatedVertices(start, end)) {
                    g_.CompressVertex(end);
                }
                g_.CompressVertex(start);
            }
        }
    }
};

//class EasyContigStream: public io::DelegatingReaderWrapper<io::SingleRead> {
//private:
//    io::Reader raw_stream_;
//    io::RCReaderWrapper<io::SingleRead> rc_stream_;
//    io::PrefixAddingReaderWrapper prefix_stream_;
//public:
//    EasyContigStream(const string& filename, const string& prefix) :
//            raw_stream_(filename), rc_stream_(raw_stream_), prefix_stream_(
//                    rc_stream_, prefix) {
//        Init(prefix_stream_);
//    }
//};

//template<class gp_t>
//class AssemblyComparer {
//private:
//    typedef typename gp_t::graph_t Graph;
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    typedef BasicSequenceMapper<Graph, typename gp_t::seq_t> Mapper; // gp_t::k_value + 1
//
//    gp_t gp_;
//    ColorHandler<Graph> coloring_;
//    io::RCReaderWrapper<io::SingleRead> rc_stream1_;
//    io::RCReaderWrapper<io::SingleRead> rc_stream2_;
//    string name1_;
//    io::PrefixAddingReaderWrapper stream1_;
//    string name2_;
//    io::PrefixAddingReaderWrapper stream2_;
//    bool untangle_;
//
////    void WriteMagicLocality() {
////        LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
////        EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);
////
////        CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
////        make_dir("/home/snurk/gingi_dbg");
////        WriteComponentsAlongPath(gp_.g, labeler, "/home/snurk/gingi_dbg/path.dot", /*split_length*/1000, /*vertex_number*/15
////                , (*MapperInstance(gp_)).MapSequence((!gp_.genome).Subseq(2024608, 2067372)), *ConstructBorderColorer(gp_.g, coloring_));
////    }
//
//    template<class gp_t2>
//    void UniversalSaveGP(
//            const gp_t2& gp/*, const visualization::graph_colorer::GraphColorer<typename gp_t2::graph_t> coloring*/,
//            const string& filename) {
//        typename PrinterTraits<Graph>::Printer printer(gp.g);
//        INFO("Saving graph to " << filename);
//        printer.saveGraph(filename);
//        printer.saveEdgeSequences(filename);
//        printer.savePositions(filename, gp.edge_pos);
////        SaveColoring(gp.g
////                , coloring
////                , filename);
//
////        LengthIdGraphLabeler<Graph> labeler(gp.g);
////        WriteSimple(gp.g, labeler, filename + ".dot");
//    }
//
//    void SaveOldGraph(const string& path) {
//        INFO("Saving graph to " << path);
//        PrintGraphPack(path, gp_);
////        LengthIdGraphLabeler<Graph> labeler(gp_.g);
////        WriteToDotFile(gp_.g, labeler, path + ".dot");
//    }
//
//    template<class gp_t2>
//    void ProduceResults(gp_t2& gp, const ColorHandler<Graph>& coloring,
//            const string& output_folder, bool detailed_output) {
////         INFO("Removing unnecessary edges");
////         DeleteVioletEdges(gp.g, coloring);
//
//// //        if (detailed_output) {
//// //            PrintColoredGraph(gp.g, coloring, gp.edge_pos,
//// //                    output_folder + "initial_pics/purple_removed.dot");
//// //            UniversalSaveGP(gp, output_folder + "saves/purple_removed");
//// //        }
//
//// //        ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
//// //        BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
////         INFO("Counting stats, outputting pictures");
////         BPGraphStatCounter<Graph> counter(gp.g, coloring, output_folder);
////         LengthIdGraphLabeler<Graph> labeler(gp.g);
////         counter.CountStats(labeler, detailed_output);
//    }
//
//    void PrepareDirs(const string& output_folder, bool detailed_output) {
//        DIR *dp;
//        if ((dp = opendir(output_folder.c_str())) == NULL) {
//            INFO("Dir " + output_folder + " did not exist, creating");
//        } else {
//            INFO("Dir " + output_folder + " purged");
//            remove_dir(output_folder);
//        }
//        utils::MakeDirPath(output_folder);
//        if (detailed_output) {
//            make_dir(output_folder + "initial_pics/");
//            make_dir(output_folder + "saves/");
//            make_dir(output_folder + "purple_edges_pics/");
//        }
//    }
//
//public:
//
//    AssemblyComparer(size_t k_value, io::IReader<io::SingleRead> &stream1,
//            io::IReader<io::SingleRead> &stream2, const string& name1,
//            const string& name2, bool untangle = false,
//            const Sequence& reference = Sequence()) :
//            gp_(k_value, "tmp", reference, 200, true), coloring_(gp_.g, 2), rc_stream1_(
//                    stream1), rc_stream2_( // TODO dir
//                    stream2), name1_(name1), stream1_(rc_stream1_, name1), name2_(
//                    name2), stream2_(rc_stream2_, name2), untangle_(untangle) {
//    }
//
//    void CompareAssemblies(const string& output_folder, bool detailed_output =
//            true, bool one_many_resolve = false,
//            const string& add_saves_path = "") {
////        VERIFY(gp_.genome.size() > 0);
//        //todo ???
//        stream1_.reset();
//        stream2_.reset();
//
//        PrepareDirs(output_folder, detailed_output);
//        vector<ContigStream*> stream_vec = { &stream1_, &stream2_ };
//        ContigStreams streams(stream_vec, false);
//
//    CoordinatesHandler<Graph> coordinates_handler;
//        ConstructColoredGraph<gp_t>(gp_, coloring_, coordinates_handler, streams);
//
//        if (gp_.genome.size() > 0) {
//            INFO("Filling ref pos " << gp_.genome.size());
////            FillPos(gp_, gp_.genome, "ref_0");
////            FillPos(gp_, !gp_.genome, "ref_1");
//
////            SimpleInDelAnalyzer<Graph> del_analyzer(gp_.g, coloring_,
////                    gp_.edge_pos,
////                    (*MapperInstance < gp_t > (gp_)).MapSequence(gp_.genome).simple_path(),
////                    kRedColorSet, output_folder);
////            del_analyzer.Analyze();
//
////            AlternatingPathsCounter<Graph> alt_count(gp_.g, coloring);
////            alt_count.CountPaths();
//
////            ContigBlockStats<Graph, Mapper> block_stats(gp_.g, gp_.edge_pos,
////                    *MapperInstance(gp_), gp_.genome, stream1_);
////            block_stats.Count();
//
////            MissingGenesAnalyser<Graph, Mapper> missed_genes(gp_.g, coloring_,
////                    gp_.edge_pos, gp_.genome, *MapperInstance(gp_),
////                    vector<pair<bool, pair<size_t, size_t>>> {
////                        make_pair(true, make_pair(260354, 260644)),
////                        make_pair(true, make_pair(300641, 300904)),
////                        make_pair(true, make_pair(300904, 301920)),
////                        make_pair(true, make_pair(301917, 302348)),
////                        make_pair(true, make_pair(260354, 260644)),
////                        make_pair(true, make_pair(300641, 300904)),
////                        make_pair(true, make_pair(300904, 301920)),
////                        make_pair(true, make_pair(301917, 302348)),
////                        make_pair(true, make_pair(302449, 304752)),
////                        make_pair(true, make_pair(263821, 264594)),
////                        make_pair(true, make_pair(265025, 265726)),
////                        make_pair(true, make_pair(265740, 266951))
////                    }
////                    , output_folder + "missed_genes/");
////
////            missed_genes.Analyze();
//        }
//
//        ////////////
////        WriteMagicLocality();
//        ////////////
//
////        2339834
////        INFO("Removing gaps");
////        GapsRemover<Graph> gaps_remover(gp_.g, coloring, kBlueColorSet, 700);
////        gaps_remover.RemoveGaps();
////        INFO("Gaps removed");
//
//        if (boost::starts_with(name1_, "idba")) {
//            IDBADiffAnalyzer<gp_t> diff_analyzer(gp_, coloring_, name1_, name2_,
//                    output_folder + "/idba_analysis/");
//            diff_analyzer.Analyze(stream1_, stream2_);
//        }
//
//        if (one_many_resolve) {
//            VERIFY(!untangle_);
//            RestrictedOneManyResolver<Graph> resolver(gp_.g, coloring_,
//                    kVioletColorSet);
//            resolver.Resolve();
//        }
//
//        if (detailed_output) {
//            if (gp_.genome.size() > 0) {
//                PrintColoredGraphAlongRef(gp_, coloring_, // gp_.edge_pos,  TODO why not corresponding???
//                        //gp_.genome,
//                        output_folder + "initial_pics/colored_split_graph.dot");
//            } else {
//                PrintColoredGraph(gp_.g, coloring_, gp_.edge_pos,
//                        output_folder + "initial_pics/colored_split_graph.dot");
//            }
//
//            if (add_saves_path != "") {
//                UniversalSaveGP(gp_, //coloring,
//                        add_saves_path);
//                SaveColoring(gp_.g, coloring_, add_saves_path);
//                //PrintColoredGraphWithColorFilter(gp_.g, coloring_, gp_.edge_pos,
//            //            add_saves_path + ".dot");
//            }
//            UniversalSaveGP(gp_, //coloring,
//                    output_folder + "saves/colored_split_graph");
//            SaveColoring(gp_.g, coloring_,
//                    output_folder + "saves/colored_split_graph");
//            //PrintColoredGraphWithColorFilter(gp_.g, coloring_, gp_.edge_pos,
//            //        output_folder + "saves/colored_split_graph.dot");
//        }
//
//        if (untangle_) {
//            VERIFY(false);
////            INFO("Untangling graph");
////            bp_graph_pack<typename gp_t::graph_t> untangled_gp(gp_t::k_value);
////            UntangledGraphConstructor<gp_t> untangler(gp_, coloring,
////                    untangled_gp, stream1_, stream2_);
////            //todo ???
////            //        SimplifyGraph(untangled_gp.g);
////
////            if (detailed_output) {
////                PrintColoredGraph(untangled_gp.g, untangled_gp.coloring,
////                        untangled_gp.edge_pos,
////                        output_folder + "initial_pics/untangled_graph.dot");
////                UniversalSaveGP(untangled_gp, //untangled_gp.coloring,
////                        output_folder + "saves/untangled_graph");
////            }
////
////            ProduceResults(untangled_gp, untangled_gp.coloring, output_folder,
////                    detailed_output);
//        } else {
////            INFO("Analyzing gaps");
////            GapComparativeAnalyzer<Graph> gap_analyzer(gp_.g, coloring,
////                    gp_.edge_pos);
////            gap_analyzer.ReportPotentialGapsCloses(
////                    output_folder + "gap_closing_edges/");
//
//            //trivial breakpoints
//            string bp_folder = output_folder + "breakpoints/";
//            make_dir(bp_folder);
//            TrivialBreakpointFinder<Graph> bp_finder(gp_.g, coloring_,
//                    gp_.edge_pos);
//            bp_finder.FindBreakPoints(bp_folder);
//
//            //possible rearrangements
//            string rearr_folder = output_folder + "rearrangements/";
//            make_dir(rearr_folder);
//            SimpleRearrangementDetector<gp_t> rearr_det(gp_, coloring_, "tdc_",
//                    rearr_folder);
//            rearr_det.Detect();
//
//            ProduceResults(gp_, coloring_, output_folder, detailed_output);
//        }
//    }
//
//private:
//    DECL_LOGGER("AssemblyComparer")
//    ;
//};

//template<size_t k, size_t K, class BuildSeq>
//void RunBPComparison(ContigStream& raw_stream1, ContigStream& raw_stream2,
//        const string& name1, const string& name2, bool refine, bool untangle,
//        const string& output_folder, bool detailed_output = true, size_t delta =
//                5, Sequence reference = Sequence(),
//        const string& add_saves_path = "") {
//    static double lol_time = 0;
//    cap::utils::add_time(lol_time, -1);
//
//    io::SplittingWrapper stream1(raw_stream1);
//    io::SplittingWrapper stream2(raw_stream2);
//
//    typedef debruijn_graph::graph_pack<
//    /*Nonc*/debruijn_graph::ConjugateDeBruijnGraph, BuildSeq> comparing_gp_t;
//
//    if (refine) {
//        typedef graph_pack<ConjugateDeBruijnGraph, BuildSeq> refining_gp_t;
//        refining_gp_t refining_gp(k, "tmp");
//        io::VectorReader<io::SingleRead> genome_stream(
//                io::SingleRead("genome", reference.str()));
//        ContigStreamsPtr streams_ptr = make_shared<ContigStreams>(vector<ContigStream*>{&stream1, &stream2, &genome_stream}, false);
//
//        ConstructGPForRefinement(refining_gp, streams_ptr, delta);
//
//        io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1,
//                GraphReadCorrectorInstance(refining_gp.g,
//                        *MapperInstance(refining_gp)));
//        io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2,
//                GraphReadCorrectorInstance(refining_gp.g,
//                        *MapperInstance(refining_gp)));
//        io::ModifyingWrapper<io::SingleRead> reference_stream(genome_stream,
//                GraphReadCorrectorInstance(refining_gp.g,
//                        *MapperInstance(refining_gp)));
//
//        reference_stream.reset();
//        AssemblyComparer<comparing_gp_t> comparer(K, refined_stream1,
//                refined_stream2, name1, name2, untangle,
//                ReadSequence(reference_stream));
//        comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/
//                false, add_saves_path);
//    } else {
//        AssemblyComparer<comparing_gp_t> comparer(K, stream1, stream2, name1,
//                name2, untangle, reference);
//        comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/
//                false, add_saves_path);
//    }
//
//    cap::utils::add_time(lol_time, +1);
//    INFO("LOL_TIME:: " << lol_time);
//}

//template<size_t k, size_t K>
//void RunBPComparison(const Sequence& ref, ContigStream& stream,
//        const string& name1, const string& name2, bool refine, bool untangle,
//        const string& output_folder, bool detailed_output = true, size_t delta =
//                5) {
//    io::VectorReader<io::SingleRead> ref_stream(
//            io::SingleRead(name1, ref.str()));
//    RunBPComparison<k, K>(ref_stream, stream, name1, name2, refine, untangle,
//            output_folder, detailed_output, delta);
//}

//template<size_t k, size_t K>
//void RunBPComparison(const Sequence& s1, const Sequence& s2,
//        const string& name1, const string& name2, bool refine, bool untangle,
//        const string& output_folder, bool detailed_output = true) {
//    io::VectorReader<io::SingleRead> stream(io::SingleRead(name2, s2.str()));
//    RunBPComparison<k, K>(s1, stream, name1, name2, refine, untangle,
//            output_folder, detailed_output);
//}
//
//template<size_t k, size_t K>
//void RunBPComparison(const Sequence& ref, const vector<Sequence>& contigs,
//        const string& name1, const string& name2, bool refine, bool untangle,
//        const string& output_folder, bool detailed_output = true) {
//    io::VectorReader<io::SingleRead> stream(MakeReads(contigs));
//    RunBPComparison<k, K>(ref, stream, name1, name2, refine, untangle,
//            output_folder, detailed_output);
//}

//template<size_t k, class BuildSeq>
//void CompareGenomes(const Sequence& genome_1, const Sequence& genome_2,
//        const string& output_dir) {
//    INFO("Genome comparison started");
//    io::VectorReader<io::SingleRead> stream1(
//            io::SingleRead("", genome_1.str()));
//    io::VectorReader<io::SingleRead> stream2(
//            io::SingleRead("", genome_2.str()));
//    typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, BuildSeq> comparing_gp_t; // k
//    INFO("Running assembly comparer");
//    AssemblyComparer<comparing_gp_t> comparer(k, stream1, stream2, "strain1",
//            "strain2", /*untangle*/false);
//    comparer.CompareAssemblies(output_dir, /*detailed_output*/true, /*on_many_resolve*/
//    false);
//    INFO("Finished");
//}

template<class gp_t>
void ThreadAssemblies(const string& base_saves, ContigStream& base_assembly,
        const string& base_prefix, ContigStream& assembly_to_thread,
        const string& to_thread_prefix, const string& output_dir) {
    typedef typename gp_t::graph_t Graph;
    gp_t gp;
//        ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, base_assembly);
    ScanGraphPack(base_saves, gp);
    base_assembly.reset();
    visualization::position_filler::FillPos(gp, base_assembly, base_prefix);
    visualization::position_filler::FillPos(gp, assembly_to_thread, to_thread_prefix);

    visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
    visualization::graph_labeler::StrGraphLabeler<Graph> str_labeler(gp.g);
    visualization::graph_labeler::CompositeLabeler<Graph> labeler(pos_labeler, str_labeler);

    auto mapper = MapperInstance(gp);

    assembly_to_thread.reset();
    io::SingleRead read;
    while (!assembly_to_thread.eof()) {
        assembly_to_thread >> read;
        make_dir(output_dir + read.name());
        WriteComponentsAlongPath(gp.g, labeler,
                output_dir + read.name() + "/.dot", /*split_edge_length*/400,
                mapper->MapSequence(read.sequence()),
                Path<typename Graph::EdgeId>(), Path<typename Graph::EdgeId>(),
                true);
    }
}

template<class gp_t>
void RunMultipleGenomesVisualization(size_t k_visualize,
        vector<pair<std::string, std::string> > genomes_paths,
        std::string output_folder) {
    typedef typename gp_t::graph_t Graph;

    utils::MakeDirPath(output_folder);

    gp_t gp(k_visualize, "tmp", 0, Sequence(), 200);
    ColorHandler<Graph> coloring(gp.g, genomes_paths.size());
  CoordinatesHandler<Graph> coordinates_handler;

    // ContigStream -> SplittingWrapper -> RCReaderWrapper -> PrefixAddingReaderWrapper

  ContigStreams streams;
  for (auto it = genomes_paths.begin(); it != genomes_paths.end(); ++it) {
    streams.push_back(make_shared<io::FileReadStream>(it->second));
  }

  ContigStreams rc_wrapped = io::RCWrap(streams);

  ConstructColoredGraph(gp, coloring, coordinates_handler, rc_wrapped);

    ofstream indel_event_logger(output_folder + "/indel_events");

//  UnversalSaveGP(gp, output_folder + "/colored_split_graph");
//  SaveColoring(gp.g, coloring, output_folder + "/colored_split_graph");
    //PrintColoredGraphWithColorFilter(gp.g, coloring, gp.edge_pos,
    //        output_folder + "/colored_split_graph.dot");
}

}

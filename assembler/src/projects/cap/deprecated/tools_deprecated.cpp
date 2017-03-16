//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


//Gingi block

//BOOST_AUTO_TEST_CASE( MaskDiffsForGingi ) {
//    MaskDifferencesAndSave(vector<string> {
//            "/home/snurk/Dropbox/gingi/jeff.fasta",
//            "/home/snurk/Dropbox/gingi/TDC60.fasta" }, vector<string> { "jeff",
//            "tdc60" }, "assembly_comp/gingi_diff_mask/", k<15>(), k<21>(),
//            k<55>()/*, k<101>(), k<201>()*/);
//}

//BOOST_AUTO_TEST_CASE( ClearGingiGenome ) {
//    Clear<201>("assembly_comp/gingi_diff_mask/tdc60.fasta",
//        "assembly_comp/gingi_diff_mask/tdc60_cl.fasta");
//}
//
//BOOST_AUTO_TEST_CASE( ClearJeffAssembly ) {
//    Clear<201>("assembly_comp/gingi_diff_mask/jeff.fasta",
//        "assembly_comp/gingi_diff_mask/jeff_cl.fasta");
//}
// BOOST_AUTO_TEST_CASE
// ( AssemblyRefComparison ) {
//     static const size_t K = 55;
//     typedef debruijn_graph::graph_pack<
//     /*Nonc*/debruijn_graph::ConjugateDeBruijnGraph, K> gp_t;
//     typedef gp_t::graph_t Graph;
//     typedef Graph::EdgeId EdgeId;
//     typedef Graph::VertexId VertexId;
//     typedef NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> Mapper;

// //    EasyContigStream stream_1("/home/snurk/Dropbox/gingi/jeff.fasta");
// //    EasyContigStream stream_2("/home/snurk/Dropbox/gingi/TDC60.fasta");
// //    string ref = "/home/snurk/Dropbox/gingi/TDC60.fasta";
// //    EasyContigStream stream_1("assembly_comp/gingi_diff_mask/jeff_cl.fasta");
// //    EasyContigStream stream_2("assembly_comp/gingi_diff_mask/tdc60_cl.fasta");
// //    string ref = "assembly_comp/gingi_diff_mask/tdc60_cl.fasta";
//     EasyContigStream stream_1("/home/snurk/Dropbox/lab/mrsa/MRSA_RCH_S60.fasta",
//             "s60_");
//     EasyContigStream stream_2(
//             "/home/snurk/Dropbox/lab/mrsa/USA300_FPR3757.fasta", "usa300_");
// //    EasyContigStream stream_1("assembly_comp/gingi_diff_mask/jeff.fasta",
// //            "jeff_");
// //    EasyContigStream stream_2("assembly_comp/gingi_diff_mask/tdc60.fasta",
// //            "tdc_");

//     string ref = "/home/snurk/Dropbox/lab/mrsa/USA300_FPR3757.fasta";
// //    string ref = "assembly_comp/gingi_diff_mask/tdc60.fasta";
//     string output_folder = "assembly_comp/s60_usa300_" + std::to_string(K) + "/";
//     remove_dir(output_folder);
//     make_dir(output_folder);

//     int br_delta = -1;
//     gp_t gp(ReadGenome(ref), 200, true);
//     ColorHandler<Graph> coloring(gp.g);

//     vector<ContigStream*> streams = { &stream_1, &stream_2 };
//     ConstructColoredGraph(gp, coloring, streams, false, br_delta);

// //    INFO("Filling ref pos " << gp.genome.size());
// //            visualization::position_filler::FillPos(gp_, gp_.genome, "ref_0");
// //            visualization::position_filler::FillPos(gp_, !gp_.genome, "ref_1");

// //Indels
// //    make_dir(output_folder + "indels/");
// //    SimpleInDelAnalyzer<Graph> del_analyzer(gp.g, coloring, gp.edge_pos,
// //            (*MapperInstance(gp)).MapSequence(gp.genome).simple_path(),
// //            edge_type::red, output_folder + "indels/");
// //    del_analyzer.Analyze();

// //Alternating paths
// //            AlternatingPathsCounter<Graph> alt_count(gp_.g, coloring);
// //            alt_count.CountPaths();

// //Block stats
// //            ContigBlockStats<Graph, Mapper> block_stats(gp_.g, gp_.edge_pos,
// //                    *MapperInstance(gp_), gp_.genome, stream1_);
// //            block_stats.Count();

// //    Missing genes
// //    MissingGenesAnalyser<Graph, Mapper> missed_genes(gp.g, coloring,
// //            gp.edge_pos, gp.genome, *MapperInstance(gp),
// //            vector<pair<bool, pair<size_t, size_t>>> {
// //            make_pair(/*true*/false, make_pair(416000, 430000)),
// //            make_pair(/*true*/false, make_pair(1513000, 1518000)),
// //            make_pair(/*true*/false, make_pair(260354, 260644)),
// //            make_pair(/*true*/false, make_pair(300641, 300904)),
// //            make_pair(/*true*/false, make_pair(300904, 301920)),
// //            make_pair(/*true*/false, make_pair(301917, 302348)),
// //            make_pair(/*true*/false, make_pair(260354, 260644)),
// //            make_pair(/*true*/false, make_pair(300641, 300904)),
// //            make_pair(/*true*/false, make_pair(300904, 301920)),
// //            make_pair(/*true*/false, make_pair(301917, 302348)),
// //            make_pair(/*true*/false, make_pair(302449, 304752)),
// //            make_pair(/*true*/false, make_pair(263821, 264594)),
// //            make_pair(/*true*/false, make_pair(265025, 265726)),
// //            make_pair(/*true*/false, make_pair(265740, 266951))
// //        }
// //        , output_folder + "missed_genes/");
// //    missed_genes.Analyze();

// //        2339834
// ////////////
// //    WriteMagicLocality();
// ////////////

// //possible rearrangements
// //        string rearr_folder = output_folder + "rearrangements/";
// //        make_dir(rearr_folder);
// //        SimpleRearrangementDetector<gp_t> rearr_det(gp_, coloring_, "tdc_",
// //                rearr_folder);
// //        rearr_det.Detect();

// //print graph
//     make_dir(output_folder + "initial_pics");
//     PrintColoredGraphAlongRef(gp, coloring,
//             output_folder + "initial_pics/colored_split_graph.dot");

//     //reference correction
//     SimpleInDelCorrector<Graph> corrector(gp.g, coloring,
//             (*MapperInstance(gp)).MapSequence(gp.genome).simple_path().sequence(), /*genome_color*/
//             kBlueColorSet, /*assembly_color*/kRedColorSet);
//     corrector.Analyze();

//     //trivial breakpoints
//     string bp_folder = output_folder + "breakpoints/";
//     make_dir(bp_folder);
//     TrivialBreakpointFinder<Graph> bp_finder(gp.g, coloring, gp.edge_pos);
//     bp_finder.FindBreakPoints(bp_folder);

//     //make saves
//     make_dir(output_folder + "saves");
//     string filename = output_folder + "saves/graph";
//     PrinterTraits<Graph>::Printer printer(gp.g);
//     INFO("Saving graph to " << filename);
//     printer.saveGraph(filename);
//     printer.saveEdgeSequences(filename);
//     printer.savePositions(filename, gp.edge_pos);
//     SaveColoring(gp.g, coloring, filename);
// }

//End of gingi block

//BOOST_AUTO_TEST_CASE( RefVSAssemblyComparison ) {
//    static const size_t k = 55;
//    static const size_t K = 55;
//    Sequence ref = ReadGenome("/home/snurk/MRSA/USA300_FPR3757.fasta");
//    io::Reader contig_stream("/home/snurk/MRSA/MRSA_RCH_I56.fasta");
//    string folder = "mrsa_comp/RCH_I56/";
//    make_dir(folder);
//    RunBPComparison<k, K>(
//            ref,
//            contig_stream,
//            "ref",
//            "assembly",
//            true/*refine*/,
//            false/*untangle*/,
//            folder,
//            true/*detailed_output*/,
//            20);
//}

//BOOST_AUTO_TEST_CASE( TwoAssemblyComparison ) {
//    static const size_t k = 19;
//    static const size_t K = 55;
////    static const size_t K = 57;
////    static const size_t K = 53;
//
////    io::Reader stream_1("/home/snurk/gingi/2.fasta");
////    io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//    io::Reader stream_1("/home/anton/idba_compare/idba.fasta");
//    io::Reader stream_2("/home/anton/idba_compare/hammer21_dis_tuned_simpl_try_improve.fasta");
//    string ref = "/home/anton/idba_compare/MG1655-K12.fasta";
//    string folder = "/home/anton/idba_compare/hammer21_dis_tuned_simpl_vs_idba/";
////    string folder = "assembly_comp/gingi_new_3_vs_jeff/";
//    make_dir(folder);
//
//    RunBPComparison<k, K>(
//        stream_1,
//        stream_2,
//        "idba_",
//        "k21ts_",
//        true/*refine*/,
//        false/*untangle*/,
//        folder,
//        true/*detailed_output*/,
//        5/*delta*/,
//        ReadGenome(ref));
//}

//BOOST_AUTO_TEST_CASE( TwoAssemblyComparison ) {
//    static const size_t k = 19;
//    static const size_t K = 55;
////    static const size_t K = 57;
////    static const size_t K = 53;
//
////    io::Reader stream_1("/home/snurk/gingi/2.fasta");
////    io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
////        io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE2_BH.fasta");
////    io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE3_BH.fasta");
//    io::Reader stream_2("/home/snurk/gingi/jeff.fasta");
//
////    io::Reader stream_2("/home/snurk/gingi/PGINGIVALIS_LANE2_BH.fasta");
//
//    io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE3_BH.fasta");
////    io::Reader stream_2("/home/snurk/gingi/lane2_evsc.fasta");
//
////    string folder = "assembly_comp/gingi_new_3_vs_new_2/";
//    string folder = "assembly_comp/gingi_new_3_vs_jeff/";
//    make_dir(folder);
//
//    RunBPComparison<k, K>(
//        stream_1,
//        stream_2,
////        "2",
////        "jeff",
//        "3_new_",
////        "2_new_",
//        "jeff_",
//        true/*refine*/,
//        false/*untangle*/,
//        folder,
//        true/*detailed_output*/);
//}

//BOOST_AUTO_TEST_CASE( AssemblyRefComparison ) {
//    static const size_t k = 21;
//    static const size_t K = 201/*55*//*201*/;
////    static const size_t K = 57;
////    static const size_t K = 53;
//
////    io::Reader stream_1("/home/snurk/gingi/2.fasta");
////    io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//    io::Reader stream_1("/home/snurk/Dropbox/gingi/jeff.fasta");
//    io::Reader stream_2("/home/snurk/Dropbox/gingi/TDC60.fasta");
//    string ref = "/home/snurk/Dropbox/gingi/TDC60.fasta";
//
////    string folder = "assembly_comp/gingi_jeff_vs_tdc60_55/";
//    string folder = "assembly_comp/gingi_jeff_vs_tdc60_201/";
//    make_dir(folder);
//
//    RunBPComparison<k, K>(
//        stream_1,
//        stream_2,
//        "jeff_",
//        "tdc_",
//        true/*refine*/,
//        false/*untangle*/,
//        folder,
//        true/*detailed_output*/,
//        5/*delta*/,
//        ReadGenome(ref));
//}

//BOOST_AUTO_TEST_CASE( IDBA_vs_SPADES ) {
//    static const size_t k = 19;
//    static const size_t K = 55;
////    static const size_t K = 57;
////    static const size_t K = 53;
//
////    io::Reader stream_1("/home/snurk/gingi/2.fasta");
////    io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//    io::Reader stream_1("/home/snurk/idba_comp/idba-contig-100.fa");
//    io::Reader stream_2("/home/snurk/idba_comp/k21nodiscard.fasta");
//    string ref = "/home/snurk/idba_comp/MG1655-K12.fasta";
//
//    string folder = "/home/snurk/idba_comp/results/";
//    make_dir(folder);
//
//    RunBPComparison<k, K>(
//        stream_1,
//        stream_2,
//        "idba_",
//        "bh21_",
//        true/*refine*/,
//        false/*untangle*/,
//        folder,
//        true/*detailed_output*/,
//        5/*delta*/,
//        ReadGenome(ref));
//}



/*
BOOST_AUTO_TEST_CASE( TwoStrainComparisonWR ) {
    INFO("Running comparison of two strains");

    make_dir("bp_graph_test");

    std::string base_dir = "/Users/valich/Dropbox/mrsa/";
    std::string genome_path1 = "/smallnas/yana/X5-l-velvet-scaff.closed.fasta",
                genome_path2 = "/smallnas/yana/X5_results/scaffolds.fasta";

    pair<Sequence, Sequence> genomes = CorrectGenomes<55>(CorrectGenomes<21>(
            ReadGenome(genome_path1),
            ReadGenome(genome_path2)), 200);

    INFO("Genomes ready");

    CompareGenomes<77>(genomes.first, genomes.second, "bp_graph_test/two_strain_comp_wr/");
    INFO("Finished");
}
*/
// inline void StrainComparisonWOR(const string& strain_1, const string& strain_2, const string& output_folder) {
//     make_dir("bp_graph_test");
//     INFO("Running comparison of two strains");
//     pair<Sequence, Sequence> genomes = CorrectGenomes<55>(TotallyClearGenomes<55>(CorrectGenomes<21>(ReadGenome(strain_1)
//             , ReadGenome(strain_2))), 30);
// //    genomes = TotallyClearGenomes<701>(genomes);
//     VERIFY(CheckNoRepeats<301>(genomes.first));
//     VERIFY(CheckNoRepeats<301>(genomes.second));
//     INFO("Genomes ready");

//     CompareGenomes<701>(genomes.first, genomes.second, output_folder);
// }

// BOOST_AUTO_TEST_CASE( TwoStrainComparisonWOR ) {
//     string base_dir = "/Users/valich/Dropbox/mrsa/";
//     StrainComparisonWOR(base_dir + "/MRSA_RCH_I56.fasta"
// , base_dir + "MRSA_RCH_S60.fasta", "bp_graph_test/two_strain_comp_wo_repeats/");
// }

//BOOST_AUTO_TEST_CASE( CompareAllMRSA ) {
//    string mrsa_root = "/home/snurk/MRSA/more_strains/";
//    ifstream stream;
//    stream.open(mrsa_root + "list.txt");
//    string s1;
//    string s2;
//    VERIFY(!stream.eof());
//    stream >> s1;
//    while (!stream.eof()) {
//        stream >> s2;
//        StrainComparisonWOR(mrsa_root + s1 + ".fasta", mrsa_root + s2 + ".fasta"
//                , mrsa_root + "results/" +  s1 + "_vs_" + s2 + "/");
//    }
//    stream.close();
//}

//
//BOOST_AUTO_TEST_CASE( TwoStrainComparisonFirstWOR ) {
//    make_dir("bp_graph_test");
//    INFO("Running comparison of two strains");
//    pair<Sequence, Sequence> genomes = CorrectGenomes<21>(ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz")
//            , ReadGenome("data/input/E.coli/DH10B-K12.fasta"));
//    genomes = CorrectGenomes<55>(genomes, 200);
//    genomes.first = ClearGenome<55>(genomes.first);
//    genomes = CorrectGenomes<55>(genomes, 200);
//    VERIFY(CheckNoRepeats<301>(genomes.first));
//    INFO("Genomes ready");
//
//    CompareGenomes<701>(genomes.first, genomes.second, "bp_graph_test/two_strain_comp_first_wo_repeats/");
//}

//BOOST_AUTO_TEST_CASE( StrainVSRepeatGraphComparison ) {
//    static const size_t repeat_clearing_k = 55;
//    static const size_t repeat_graph_k = 101;
//    static const size_t refining_k1 = 25;
//    static const size_t refining_k2 = 151;
//    static const size_t bp_k = 301;
//
//    make_dir("bp_graph_test");
//    INFO("Running comparison of strain vs repeat graph for other strain");
//    Sequence genome1 = ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz");
//    Sequence genome2 = ReadGenome("data/input/E.coli/DH10B-K12.fasta");
//
//    typedef graph_pack<ConjugateDeBruijnGraph, repeat_graph_k> repeat_gp_t;
//    vector<Sequence> repeat_graph_edges = RepeatGraphEdges<repeat_gp_t>(genome2);
//
//    pair<Sequence, vector<Sequence>> refined_data1 = RefineData<refining_k2>(RefineData<refining_k1>(
//            make_pair(genome1, repeat_graph_edges)));
//
//    pair<Sequence, vector<Sequence>> cleared = Clear<repeat_clearing_k>(refined_data1.first, refined_data1.second);
//
//    pair<Sequence, vector<Sequence>> refined_data2 = RefineData<bp_k>(RefineData<refining_k2>(RefineData<refining_k1>(cleared)));
//
//    io::VectorReader<io::SingleRead> final_stream1(
//            io::SingleRead("first", refined_data2.first.str()));
//    io::VectorReader<io::SingleRead> final_stream2(MakeReads(refined_data2.second));
//
//    typedef graph_pack<ConjugateDeBruijnGraph, bp_k> comparing_gp_t;
//    INFO("Running assembly comparer");
//    AssemblyComparer<comparing_gp_t> comparer(final_stream1, final_stream2, "strain1", "strain2", /*untangle*/false);
//    comparer.CompareAssemblies("bp_graph_test/strain_vs_repeat_graph_comp/", /*detailed_output*/true);
//    INFO("Finished");
//}
//
//BOOST_AUTO_TEST_CASE( StrainVSRepeatGraphComparison2 ) {
//    static const size_t repeat_clearing_k = 55;
//    static const size_t repeat_graph_k = 201;
//    static const size_t refining_k1 = 25;
//    static const size_t refining_k2 = 151;
//    static const size_t bp_k = 201;
//
//    make_dir("bp_graph_test");
//    INFO("Running comparison of strain vs repeat graph for other strain");
//    Sequence genome1 = ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz");
//    Sequence genome2 = ReadGenome("data/input/E.coli/DH10B-K12.fasta");
//
//    typedef graph_pack<ConjugateDeBruijnGraph, repeat_graph_k> repeat_gp_t;
//    vector<Sequence> repeat_graph_edges = RepeatGraphEdges<repeat_gp_t>(genome2);
//
//    pair<Sequence, vector<Sequence>> refined_data1 = RefineData<refining_k2>(RefineData<refining_k1>(
//            make_pair(genome1, repeat_graph_edges)));
//
////    pair<Sequence, vector<Sequence>> cleared = Clear<repeat_clearing_k>(refined_data1.first, refined_data1.second);
//    Sequence cleared = ClearGenome<repeat_clearing_k>(refined_data1.first);
//
//    pair<Sequence, vector<Sequence>> refined_data2 = RefineData<bp_k>(RefineData<refining_k2>(RefineData<refining_k1>(make_pair(cleared, refined_data1.second))));
//
//    io::VectorReader<io::SingleRead> final_stream1(
//            io::SingleRead("first", refined_data2.first.str()));
//    io::VectorReader<io::SingleRead> final_stream2(MakeReads(refined_data2.second));
//
//    typedef graph_pack<ConjugateDeBruijnGraph, bp_k> comparing_gp_t;
//    INFO("Running assembly comparer");
//    AssemblyComparer<comparing_gp_t> comparer(final_stream1, final_stream2, "strain1", "strain2", /*untangle*/false);
//    comparer.CompareAssemblies("bp_graph_test/strain_vs_repeat_graph_comp2/", /*detailed_output*/true);
//    INFO("Finished");
//}

//BOOST_AUTO_TEST_CASE( ThreadingContigsOverGraph ) {
//    typedef graph_pack<ConjugateDeBruijnGraph, 55> gp_t;
//    io::EasyReader base_contigs("/home/anton/gitrep/algorithmic-biology/assembler/data/tmp/andrew_nurk.fasta");
//    io::EasyReader other_contigs("/home/anton/gitrep/algorithmic-biology/assembler/data/tmp/velvet-sc.fasta");
//    string base_saves = "/home/anton/gitrep/algorithmic-biology/assembler/data/debruijn/LBOUILLONII_QUAKE/saves/simplified_graph";
//    string output_dir = "bul_comparison/ac";
//    make_dir(output_dir);
//    ThreadAssemblies<gp_t>(base_saves, base_contigs, "spades", other_contigs, "velvet", output_dir);
//}

//BOOST_AUTO_TEST_CASE( GenerateGraphFragment ) {
//    std::string input_path = "./data/debruijn/HMP_LANE_3_0/K55/latest/saves/simplified_graph";
//    std::string output_path = "./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path";
//    size_t split_threshold = 230;
//    int int_edge_id = 5573881;
//    graph_pack<ConjugateDeBruijnGraph, 55> gp;
//    ScanGraphPack(input_path, gp);
//    //prints only basic graph structure
//    PrintGraphComponentContainingEdge(output_path, gp.g,
//            split_threshold, int_edge_id);
//
//    //long way to write to dot file
//    Graph g(55);
//    ScanBasicGraph(output_path, g);
//    total_labeler_graph_struct graph_struct(g, (const EdgesPositionHandler<Graph>*)0);
//    total_labeler tot_lab(&graph_struct);
//    WriteToDotFile(g,
//            tot_lab, output_path + ".dot",
//            "mygraph", Path<EdgeId>(), Path<EdgeId>());
//}

/*
BOOST_AUTO_TEST_CASE( GapComparativeAnalysis ) {
    std::string strain1 = "/smallnas/yana/X5-l-velvet-scaff.closed.fasta",
                strain2 = "/smallnas/yana/X5_results/scaffolds_fcb.fasta";
}
*/

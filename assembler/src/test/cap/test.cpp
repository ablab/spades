//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "compare_standard.hpp"
#include "synthetic_tests.hpp"
#include "utils/logger/log_writers.hpp"
#include "pipeline/graphio.hpp"
#include <boost/test/unit_test.hpp>

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"
#include "assembly_compare.hpp"
#include "test_utils.hpp"
#include "junk_cropping_reader.hpp"

::boost::unit_test::test_suite* init_unit_test_suite(int, char*[]) {
    //logging::create_logger("", logging::L_DEBUG);
    //logging::__logger()->add_writer(make_shared<logging::console_writer>());
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);

    using namespace ::boost::unit_test;
    char module_name[] = "cap_test";

    assign_op(framework::master_test_suite().p_name.value,
              basic_cstring<char>(module_name), 0);

    return 0;
}

#include "lseq_test.hpp"

namespace cap {

inline void CheckDiffs(size_t k, const string& actual_prefix,
                       const string& etalon_prefix, const string& work_dir,
                       bool exact_match = true) {
    if (exact_match) {
        INFO("Checking differences for graphs: exact match");
        vector<string> suffixes = { "grp", "clr", "sqn", "blk" };
        for (auto suff: suffixes) {
            BOOST_CHECK_MESSAGE(
                    CheckFileDiff(actual_prefix + "." + suff, etalon_prefix + "." + suff),
                    "Check for suffix " + suff + " failed");
        }
    } else {
        INFO("Checking differences for graphs: graph isomorphism");
        //currently not testing coordinates handler at all
        ColoredGraphIsomorphismChecker<conj_graph_pack> checker(k, work_dir);
        BOOST_CHECK_MESSAGE(checker.Check(actual_prefix, etalon_prefix),
                            "GRAPHS DIFFER");
    }
}

inline void RegenerateEtalon(size_t k, const string& filename,
                             const string& etalon_root,
                             const string& work_dir) {
    SyntheticTestsRunner<RtSeq> test_runner(filename, k, etalon_root,
                                                       work_dir);
    test_runner.Run();
}

template<class Seq>
void RunTests(size_t k, const string& filename, const string& output_dir,
              const string& etalon_root, const string& work_dir,
              bool exact_match = true, const vector<size_t>& /*example_ids*/ =
                      vector<size_t>()) {
    SyntheticTestsRunner<Seq> test_runner(filename, k, output_dir, work_dir);
    vector<size_t> launched = test_runner.Run();
    for (size_t id: launched) {
        CheckDiffs(k, output_dir + std::to_string(id), etalon_root + std::to_string(id),
                   work_dir, exact_match);
    }
}

BOOST_AUTO_TEST_CASE( RegenerateEtalonTest ) {
    return;
    utils::TmpFolderFixture _("tmp");
    string input_dir = "./src/test/cap/tests/synthetic/";
    string etalon_dir = input_dir + "etalon/";
    remove_dir(etalon_dir);
    make_dir(etalon_dir);
    RegenerateEtalon(25, input_dir + "tests.xml", etalon_dir, "tmp");
}

BOOST_AUTO_TEST_CASE( SyntheticExamplesTestsRtSeq ) {
//    return;
    utils::TmpFolderFixture _("tmp");
    string input_dir = "./src/test/cap/tests/synthetic/";
    RunTests<RtSeq>(25, input_dir + "tests.xml", "tmp/",
                               input_dir + "etalon/", "tmp", /*true*/false);
}

BOOST_AUTO_TEST_CASE( SyntheticExamplesTestsLSeq ) {
//    return;
    utils::TmpFolderFixture _("tmp");
    string input_dir = "./src/test/cap/tests/synthetic/";
    RunTests<LSeq>(25, input_dir + "tests.xml", "tmp/", input_dir + "etalon/",
                   "tmp", false);
}

/*
 BOOST_AUTO_TEST_CASE( SyntheticExamplesWithErrorsTests ) {
 utils::TmpFolderFixture _("tmp");
 make_dir("bp_graph_test");
 LoadAndRunBPG<15, 25, RtSeq>("./src/test/cap/tests/synthetic_with_err/tests2.xml",
 "bp_graph_test/simulated_common_err/", "./src/test/cap/tests/synthetic_with_err/etalon/", "1_err", false);
 LoadAndRunBPG<15, 25, LSeq>("./src/test/cap/tests/synthetic_with_err_lseq/tests2.xml",
 "bp_graph_test/simulated_common_err/", "./src/test/cap/tests/synthetic_with_err_lseq/etalon/", "1_err_lseq", false);
 remove_dir("bp_graph_test");
 }
 */

BOOST_AUTO_TEST_CASE( RepeatCroppingReaderTest ) {
    ContigStreamPtr raw_reader = make_shared<io::VectorReadStream<io::SingleRead>>(MakeReads(vector<string>{
        "ACGTCacgtcTTGCA"}));
    io::SingleRead read;
    (*raw_reader) >> read;
    BOOST_CHECK_EQUAL("ACGTCacgtcTTGCA", read.GetSequenceString());
    JunkCroppingWrapper reader(raw_reader);
    reader.reset();
    reader >> read;
    BOOST_CHECK_EQUAL("ACGTCTTGCA", read.sequence().str());
    vector<pair<size_t, size_t>> etalon_ladder = {{0, 0}, {5, 5}, {5, 10}, {10, 15}};
    BOOST_CHECK_EQUAL(reader.coordinates_ladder(), etalon_ladder);
}

BOOST_AUTO_TEST_CASE( RepeatCroppingReaderTest2 ) {
    ContigStreamPtr raw_reader = make_shared<io::VectorReadStream<io::SingleRead>>(MakeReads(vector<string>{
        "acgtcACGTCNNNNNTTGCADMYNY"}));
    io::SingleRead read;
    (*raw_reader) >> read;
    BOOST_CHECK_EQUAL("acgtcACGTCNNNNNTTGCADMYNY", read.GetSequenceString());
    JunkCroppingWrapper reader(raw_reader);
    reader.reset();
    reader >> read;
    BOOST_CHECK_EQUAL("ACGTCTTGCA", read.sequence().str());
    vector<pair<size_t, size_t>> etalon_ladder = {{0, 0}, {0, 5}, {5, 10}, {5, 15}, {10, 20}, {10, 25}};
    BOOST_CHECK_EQUAL(reader.coordinates_ladder(), etalon_ladder);

    CoordinatesHandler<cap::Graph> coords;
    coords.StoreGenomeThreadManual(0, reader.coordinates_ladder());
    for (size_t i = 0; i < etalon_ladder.size(); ++i) {
        const auto &p = etalon_ladder[i];
        if (i > 0 && p.first == etalon_ladder[i - 1].first)
            continue;

        size_t orig_pos = coords.GetOriginalPos(0, coords.PreprocessCoordinates(p.first));
        size_t etalon_pos = coords.PreprocessCoordinates(p.second);
        DEBUG("get " << debug::PrintComplexPosition(orig_pos) << " etalon " << debug::PrintComplexPosition(etalon_pos));
        BOOST_CHECK_EQUAL(orig_pos, etalon_pos);
    }
}

}

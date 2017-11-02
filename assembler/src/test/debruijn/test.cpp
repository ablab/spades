//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"

#include "test_utils.hpp"

//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"
#include "path_extend_test.hpp"
#include "overlap_removal_test.hpp"
#include "overlap_analysis_test.hpp"
//#include "detail_coverage_test.hpp"
#include "histogram_test.hpp"
#include "paired_info_test.hpp"
#include "io_test.hpp"
#include "graph_alignment_test.hpp"

#define BOOST_TEST_SOURCE
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/results_collector.ipp>
#include <boost/test/impl/unit_test_log.ipp>
#include <boost/test/impl/framework.ipp>
#include <boost/test/impl/progress_monitor.ipp>
#include <boost/test/impl/execution_monitor.ipp>
#include <boost/test/impl/unit_test_parameters.ipp>
#include <boost/test/impl/unit_test_monitor.ipp>
#include <boost/test/impl/xml_log_formatter.ipp>
#include <boost/test/impl/xml_report_formatter.ipp>
#include <boost/test/impl/plain_report_formatter.ipp>
#include <boost/test/impl/junit_log_formatter.ipp>
#include <boost/test/impl/debug.ipp>
#include <boost/test/impl/test_tree.ipp>
#include <boost/test/impl/test_tools.ipp>
#include <boost/test/impl/compiler_log_formatter.ipp>
#include <boost/test/impl/results_reporter.ipp>
#include <boost/test/impl/decorator.ipp>

::boost::unit_test::test_suite*    init_unit_test_suite( int, char* [] )
{
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);

    // Fix number of threads according to OMP capabilities.
    int max_threads = std::min(4, omp_get_max_threads());
    // Inform OpenMP runtime about this :)
    omp_set_num_threads(max_threads);
    
    using namespace ::boost::unit_test;
    char module_name [] = "debruijn_test";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

    return 0;
}

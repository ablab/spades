//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "standard_base.hpp"
#include "logger/log_writers.hpp"

#include "graphio.hpp"
#include "test_utils.hpp"

//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"
#include "order_and_law_test.hpp"
#include "path_extend_test.hpp"
#include "overlap_analysis_test.hpp"
//#include "detail_coverage_test.hpp"
//fixme why is it disabled
//#include "pair_info_test.hpp"

::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
{
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);

    using namespace ::boost::unit_test;
    char module_name [] = "debruijn_test";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

    return 0;
}

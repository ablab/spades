//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"
#include "pipeline/graphio.hpp"
#include <boost/test/unit_test.hpp>

::boost::unit_test::test_suite*    init_unit_test_suite( int, char* [] )
{
    logging::logger *log = logging::create_logger("", logging::L_DEBUG);
    log->add_writer(make_shared<logging::console_writer>());
    logging:attach_logger(log);

    using namespace ::boost::unit_test;
    char module_name [] = "debruijn_tools";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

    return 0;
}

namespace debruijn_graph {

}

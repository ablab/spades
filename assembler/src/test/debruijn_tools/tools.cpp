//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "logger/log_writers.hpp"
#include "graphio.hpp"
#include <boost/test/unit_test.hpp>
#include "io/wrapper_collection.hpp"
#include "io/osequencestream.hpp"

::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
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

BOOST_AUTO_TEST_CASE( FixReference ) {
    string fn = "";
    string out_fn = "";
    auto reader = make_shared<io::NonNuclCollapsingWrapper>(make_shared<io::FileReadStream>(fn));
    io::SingleRead read;
    std::stringstream ss;
    while (!reader->eof()) {
        (*reader) >> read;
        ss << read.GetSequenceString();
    }
    io::SingleRead concat("concat", ss.str());
    io::osequencestream out(out_fn);
    out << concat;
}

}

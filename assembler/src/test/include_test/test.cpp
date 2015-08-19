//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//#define BOOST_TEST_MODULE include_test
#include "standard_base.hpp"

#include "logger/log_writers.hpp"

#include "seq_test.hpp"
#include "rtseq_test.hpp"
#include "sequence_test.hpp"
#include "quality_test.hpp"
#include "nucl_test.hpp"

::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
{
    logging::logger *log = logging::create_logger("", logging::L_DEBUG);
    log->add_writer(std::make_shared<logging::console_writer>());
    attach_logger(log);

    using namespace ::boost::unit_test;
    char module_name [] = "include_test";
    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

    return 0;
}


//todo add more tests
////#include "ireadstream_test.hpp"
//#include "online_graph_visualizer_test.hpp"
//#include "similar_test.hpp"
//#include "single_read_test.hpp"
//#include "paired_read_test.hpp"
//#include "parser_test.hpp"
//#include "fasta_fastq_gz_parser_test.hpp"
//#include "sam_bam_parser_test.hpp"
//#include "sff_parser_test.hpp"
//#include "reader_singleread_test.hpp"
//#include "reader_pairedread_test.hpp"
//#include "multifile_reader_test.hpp"
//#include "cutting_reader_wrapper_test.hpp"
//#include "rc_reader_wrapper_test.hpp"
//#include "converting_reader_wrapper_test.hpp"

//using namespace std;


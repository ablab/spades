#include "standard_base.hpp"
#include "logger/log_writers.hpp"

#include "graphio.hpp"
#include "test_utils.hpp"

//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"
#include "order_and_law_test.hpp"
//#include "pair_info_test.hpp"


::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
{
	logging::create_logger();
	logging::__logger()->add_writer(make_shared<logging::console_writer>());

    using namespace ::boost::unit_test;
	char module_name [] = "debruijn_test";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

	return 0;
}


namespace debruijn_graph {

}


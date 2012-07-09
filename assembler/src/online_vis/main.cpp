// just to check that headers from include and debruijn folders are correctly included
#include "standard.hpp"
#include "logger/log_writers.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "config_struct.hpp"
#include "io/easy_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "launch.hpp"
#include "simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "omni/distance_estimation.hpp"
#include "memory_limit.hpp"
#include "boost/archive/tmpdir.hpp"
#include "read_converter.hpp"

#include "online_pictures.hpp"


int main(int argc, char** argv) {
	BOOST_STATIC_ASSERT(debruijn_graph::K % 2 != 0);
    
    const size_t GB = 1 << 30;

    try {
        
        using namespace online_visualization;
        string cfg_filename = argv[1];
        checkFileExistenceFATAL(cfg_filename);
        
        cfg::create_instance(cfg_filename);
        std::cout << "Hello user!" << std::endl;
        limit_memory(cfg::get().max_memory * GB);
        OnlineVisualizer online_vis;
        online_vis.run();
    }
    catch (std::exception const& e)
    {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    }
    catch (...)
    {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }
    return 0;
}

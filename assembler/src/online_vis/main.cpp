// just to check that headers from include and debruijn folders are correctly included
#include "standard.hpp"
#include "logger/log_writers.hpp"

#undef INFO
#define INFO(message)                       LOG_MSG(logging::L_DEBUG , message)

#include "standard_vis.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "config_struct.hpp"
#include "io/easy_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "memory_limit.hpp"
#include "boost/archive/tmpdir.hpp"
#include "read_converter.hpp"

#include "online_pictures.hpp"

void create_console_logger(fs::path cfg_filename)
{
	using namespace logging;

	fs::path log_props_file (cfg::get().log_filename);
	if (!exists(log_props_file))
		log_props_file = fs::path(cfg_filename).parent_path() / cfg::get().log_filename;

	create_logger(exists(log_props_file) ? log_props_file.string() : "");
	__logger()->add_writer(make_shared<console_writer>());
}



int main(int argc, char** argv) {
    const size_t GB = 1 << 30;

    try {
        
        using namespace online_visualization;
        string cfg_filename = argv[1];
        checkFileExistenceFATAL(cfg_filename);
        
        cfg::create_instance(cfg_filename);

        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);
    
        create_console_logger(cfg_filename);
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

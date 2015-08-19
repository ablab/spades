// just to check that headers from include and debruijn folders are correctly included
#include "standard.hpp"
#include "cap_kmer_index.hpp"
#include "cap_logger.hpp"

#include "../online_vis/standard_vis.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "config_struct.hpp"
#include "io/io_helper.hpp"
#include "simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "memory_limit.hpp"
#include "read_converter.hpp"

#include "cap_online_visualizer.hpp"

void create_console_logger(string const& cfg_filename) {
	using namespace logging;

    string log_props_file = cfg::get().log_filename;

    if (!path::FileExists(log_props_file))
        log_props_file = path::append_path(path::parent_path(cfg_filename), cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());

    attach_logger(lg);
}

int main(int argc, char** argv) {
    const size_t GB = 1 << 30;
    try {
        using namespace online_visualization;

        VERIFY(argc >= 2);
        string cfg_filename = argv[1];
        string cap_cfg_filename = argv[2];
        path::CheckFileExistenceFATAL(cfg_filename);

        cfg::create_instance(cfg_filename);
        cap_cfg::create_instance(cap_cfg_filename);

        create_console_logger(cfg_filename);
        cout << "\ncapGAF (Graph Analysis Framework) started" << endl;
        cout << "Print help to see help (killer feature)" << endl;
        limit_memory(cfg::get().max_memory * GB);

        CapOnlineVisualizer online_vis;
        online_vis.init();
        string batch = "";
        if (argc > 3) {
            batch = string(argv[3]);
        }
        online_vis.run(batch);
    }
    catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    }
    catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }
    return 0;
}

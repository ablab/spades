#include <string>
#include "tslr_launch.hpp"
#include <projects/online_vis/vis_logger.hpp>

void load_config(const vector<string>& cfg_fns) {
    for (const auto& s : cfg_fns) {
        path::CheckFileExistenceFATAL(s);
    }

    cfg::create_instance(cfg_fns);

    if (!cfg::get().project_name.empty()) {
        make_dir(cfg::get().output_base + cfg::get().project_name);
    }

    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().tmp_dir);

    if (cfg::get().developer_mode)
        make_dir(cfg::get().output_saves);

    make_dir(cfg::get().temp_bin_reads_path);
}

void create_console_logger(string cfg_filename) {
    using namespace logging;

    string log_props_file = cfg::get().log_filename;

    if (!path::FileExists(log_props_file))
        log_props_file = path::append_path(path::parent_path(cfg_filename), cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main (int argc, char** argv) {

	try {
        string cfg_dir = path::parent_path(argv[1]);

        vector<string> cfg_fns;
        for (int i = 1; i < argc; ++i) {
            cfg_fns.push_back(argv[i]);
        }

        load_config(cfg_fns);

        create_console_logger(cfg_dir);

        for (const auto& cfg_fn : cfg_fns)
        INFO("Loading config from " << cfg_fn);


        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);

		spades::run_tslr_resolver();
	} catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
	}
	return 0;
}

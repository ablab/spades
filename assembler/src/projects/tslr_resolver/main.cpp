#include <string>
#include "tslr_launch.hpp"
#include <modules/dev_support/copy_file.hpp>
#include <projects/online_vis/vis_logger.hpp>

void load_config(string cfg_filename) {
    path::CheckFileExistenceFATAL(cfg_filename);

    cfg::create_instance(cfg_filename);

    if (!cfg::get().project_name.empty()) {
        make_dir(cfg::get().output_base + cfg::get().project_name);
    }

    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().tmp_dir);

    make_dir(cfg::get().output_dir);
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

int main (int /*argc*/, char** argv) {

	try {
		string cfg_filename = argv[1];
        load_config(cfg_filename);
		create_console_logger(cfg_filename);
        std::cout << cfg_filename << std::endl;

        string path_to_tslr_dataset = argv[2];
        std::cout << path_to_tslr_dataset << std::endl;

        string path_to_reference = argv[3];
        std::cout << path_to_reference << std::endl;
		spades::run_tslr_resolver(path_to_tslr_dataset, path_to_reference);
	} catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
	}
	return 0;
}

/*
 * Assembler Main
 */

#include "standard.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "config_struct.hpp"
#include "io/easy_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "launch.hpp"
#include "logging.hpp"
#include "simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "omni/distance_estimation.hpp"
#include "memory_limit.hpp"

#include "assembly_compare.hpp"
#include "perfcounter.hpp"

DECL_PROJECT_LOGGER("d")

void link_output(std::string const& link_name)
{
    std::string link = cfg::get().output_root + link_name;
    unlink(link.c_str());
    if (symlink(cfg::get().output_suffix.c_str(), link.c_str()) != 0)
        WARN( "Symlink to \"" << link << "\" launch failed");
}

void link_previous_run(std::string const& previous_link_name, std::string const& link_name){
   char buf[255];

   std::string link = cfg::get().output_dir + previous_link_name;
   unlink(link.c_str());
   int count = readlink((cfg::get().output_root + link_name).c_str(), buf, sizeof(buf) - 1);
   if (count >= 0){
       buf[count] = '\0';
       std::string previous_run("../");
       previous_run = previous_run + buf;
       if (symlink(previous_run.c_str(), link.c_str()) != 0) {
           DEBUG( "Symlink to \"" << link << "\" launch failed : " << previous_run);
       }
   } else {
	   DEBUG( "Symlink to \"" << link << "\" launch failed");
   }
}

struct on_exit_output_linker
{
    on_exit_output_linker(std::string const& link_name) :
            link_name_(link_name)
    {
    }

    ~on_exit_output_linker()
    {
        link_previous_run("previous", link_name_);
        link_output(link_name_);
    }

private:
    std::string link_name_;
};

void copy_configs(fs::path cfg_filename, fs::path to)
{
    using namespace debruijn_graph;

    make_dir(to);
    copy_files_by_ext(cfg_filename.parent_path(), to, ".info", true);
}

void load_config(string cfg_filename)
{
    checkFileExistenceFATAL(cfg_filename);

    fs::path tmp_folder = fs::path(tmpnam (NULL)).parent_path() / debruijn_graph::MakeLaunchTimeDirName() / ("K" + lexical_cast<string>(debruijn_graph::K));
    copy_configs(cfg_filename, tmp_folder);

    cfg_filename = (tmp_folder / fs::path(cfg_filename).filename()).string();
    cfg::create_instance(cfg_filename);

    make_dir(cfg::get().output_root);
    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().output_saves);

    fs::path path_to_copy = fs::path(cfg::get().output_dir) / "configs";
    copy_configs(cfg_filename, path_to_copy);
}

void save_info_file() {
	ofstream file(cfg::get().output_dir + "result.info");
	file << "first\t" << cfg::get().input_dir + cfg::get().ds.first << endl;
	file << "second\t" << cfg::get().input_dir + cfg::get().ds.second << endl;
	file << "reference\t" << cfg::get().ds.reference_genome_filename << endl;
	file << "contigs\t" << cfg::get().final_contigs_file << endl;
	file.close();
}

int main(int argc, char** argv)
{
	perf_counter pc;


    const size_t GB = 1 << 30;
    limit_memory(120 * GB);

    segfault_handler sh(bind(link_output, "latest"));

    try
    {
        using namespace debruijn_graph;
        load_config(argv[1]);

        on_exit_output_linker try_linker("latest");

        // check config_struct.hpp parameters
        if (K % 2 == 0)
            VERIFY_MSG(false, "K in config.hpp must be odd!\n");

        // read configuration file (dataset path etc.)
        string dataset = cfg::get().dataset_name;

        // typedefs :)
//        typedef io::EasyReader ReadStream;
//        typedef io::PairedEasyReader PairedReadStream;

        // assemble it!
        INFO("Assembling " << dataset << " dataset with K=" << debruijn_graph::K);

        debruijn_graph::assemble_genome();

        link_output("latest_success");

        save_info_file();

        INFO("Assembling " << dataset << " dataset with K=" << debruijn_graph::K << " finished");

    }
    catch (std::exception const& e)
    {
        std::cerr << "Exception caught " << e.what() << std::endl;
        print_stacktrace();
        return EINTR;
    }
    catch (...)
    {
        std::cerr << "Unknown exception caught " << std::endl;
        print_stacktrace();
        return EINTR;
    }

	INFO("Assembling time: " << pc.time_ms());

    // OK
    return 0;
}


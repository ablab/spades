/*
 * Assembler Main
 */

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
//#include <distance_estimation.hpp>

#include "memory_limit.hpp"

DECL_PROJECT_LOGGER("d")

void link_output(std::string const& link_name)
{
	string link = cfg::get().output_root + link_name;
	unlink(link.c_str());

	if (symlink(cfg::get().output_suffix.c_str(), link.c_str()) != 0)
	    WARN( "Symlink to \"" << link << "\" launch failed");
}

struct on_exit_ouput_linker
{
	on_exit_ouput_linker(std::string const& link_name)
		: link_name_(link_name)
	{
	}

	~on_exit_ouput_linker()
	{
		link_output(link_name_);
	}

private:
	std::string link_name_;
};

void print_trace()
{
	std::cout << "=== Stack Trace ===" << std::endl;

	const size_t max_stack_size = 1000;

	void* stack_pointers[max_stack_size];
	int count = backtrace(stack_pointers, max_stack_size);

	char** func_names = backtrace_symbols(stack_pointers, count);

	// Print the stack trace
	for(int i = 0; i < count; ++i)
		std::cout << func_names[i] << std::endl;

	// Free the string pointers
	free(func_names);
}

void segfault_handler(int signum)
{
	if (signum == SIGSEGV)
	{
		std::cout << "The program was terminated by segmentation fault" << std::endl;
		print_trace();

		link_output("latest");
	}

	signal(signum, SIG_DFL);
	kill  (getpid(), signum);
}

int main() {
    const size_t GB = 1 << 30;
    limit_memory(120 * GB);

	signal(SIGSEGV, segfault_handler);

    try
    {
		using namespace debruijn_graph;

		checkFileExistenceFATAL(cfg_filename);
		cfg::create_instance(cfg_filename);

	    on_exit_ouput_linker try_linker("latest");

		// check config_struct.hpp parameters
		if (K % 2 == 0)
			VERIFY_MSG(false, "K in config.hpp must be odd!\n");

		// read configuration file (dataset path etc.)
		string input_dir = cfg::get().input_dir;
		string dataset   = cfg::get().dataset_name;

		make_dir(cfg::get().output_root );
		make_dir(cfg::get().output_dir  );
		make_dir(cfg::get().output_saves);

		// typedefs :)
		typedef io::EasyReader<io::SingleRead> ReadStream;
		typedef io::EasyReader<io::PairedRead> PairedReadStream;
//		typedef io::RCReaderWrapper<io::PairedRead> RCStream;
//		typedef io::CarefulFilteringReaderWrapper<io::PairedRead> CarefulFilteringStream;

		// assemble it!
		INFO("Assembling " << dataset << " dataset");
		INFO("K = " << debruijn_graph::K);

		debruijn_graph::assemble_genome();

		on_exit_ouput_linker("latest_success");

		INFO("Assembling " << dataset << " dataset finished");
    }
    catch(std::exception const& e)
    {
    	std::cout << "Exception caught " << e.what() << std::endl;
    }
    catch(...)
    {
    	std::cout << "Unknown exception caught " << std::endl;
    }

	// OK
	return 0;
}


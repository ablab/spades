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

    on_exit_ouput_linker try_linker("latest");

	signal(SIGSEGV, segfault_handler);

    try
    {
//    	std::cout << "false = " << debruijn_graph::CheckContains<3>(Seq<3>("ACG"), Sequence("AAAAAAA")) << std::endl;
//    	std::cout << "true = " << debruijn_graph::CheckContains<3>(Seq<3>("ACG"), Sequence("AAACGAAA")) << std::endl;
//    	std::cout << "false = " << debruijn_graph::CheckContainsSubKmer<3>(Sequence("ACGTACGT"), Sequence("AAAAAAA")) << std::endl;
//    	std::cout << "true = " << debruijn_graph::CheckContainsSubKmer<3>(Sequence("TACGACGT"), Sequence("AAACGAAA")) << std::endl;
//    	throw (std::exception());

		using namespace debruijn_graph;

		checkFileExistenceFATAL(cfg_filename);
		cfg::create_instance(cfg_filename);

		// check config_struct.hpp parameters
		if (K % 2 == 0)
			VERIFY_MSG(false, "K in config.hpp must be odd!\n");

		// read configuration file (dataset path etc.)
		string input_dir = cfg::get().input_dir;
		string dataset   = cfg::get().dataset_name;

		make_dir(cfg::get().output_root );
		make_dir(cfg::get().output_dir  );
		make_dir(cfg::get().output_saves);

		string reads_filename_1 = input_dir + cfg::get().ds.first;
		string reads_filename2 = input_dir + cfg::get().ds.second;

		checkFileExistenceFATAL(reads_filename_1);
		checkFileExistenceFATAL(reads_filename2);

		// typedefs :)
		typedef io::EasyReader<io::SingleRead> ReadStream;
		typedef io::EasyReader<io::PairedRead> PairedReadStream;
//		typedef io::RCReaderWrapper<io::PairedRead> RCStream;
//		typedef io::CarefulFilteringReaderWrapper<io::PairedRead> CarefulFilteringStream;

		// read data ('reads')

		PairedReadStream /*raw_*/paired_stream(std::make_pair(reads_filename_1,reads_filename2), cfg::get().ds.IS);

//		CarefulFilteringStream filtered_stream(raw_paired_stream);
//		RCStream paired_stream(filtered_stream);


		// read data ('genome')
		string genome_filename = cfg::get().reference_genome;
		std::string genome;
		if (genome_filename.length() > 0) {
			genome_filename = input_dir + genome_filename;
			checkFileExistenceFATAL(genome_filename);
			ReadStream genome_stream(genome_filename);
			io::SingleRead full_genome;
			genome_stream >> full_genome;
			genome = full_genome.GetSequenceString().substr(0, cfg::get().ds.LEN); // cropped
		}

		// assemble it!
		INFO("Assembling " << dataset << " dataset");
		INFO("K = " << debruijn_graph::K);

		debruijn_graph::assemble_genome(paired_stream, Sequence(genome)/*, work_tmp_dir, reads*/);

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


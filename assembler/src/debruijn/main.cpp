/*
 * Assembler Main
 */
#include "launch.hpp"
#include "config.hpp"
#include "common/logging.hpp"
#include "common/simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>

namespace {

std::string MakeLaunchTimeDirName() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d.%m_%H_%M", timeinfo);
	return string(buffer);
}
}
DECL_PROJECT_LOGGER("d")

int main()
{
	// check config.hpp parameters
	if (K % 2 == 0) {
		FATAL("K in config.hpp must be odd!\n");
	}
	checkFileExistenceFATAL(CONFIG_FILENAME);

	// read configuration file (dataset path etc.)
	string input_dir = CONFIG.read<string> ("input_dir");
	string output_dir = CONFIG.read<string> ("output_dir")
			+ MakeLaunchTimeDirName() + "/";
//	std::cout << "here " << mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH| S_IWOTH) << std::endl;
	mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH| S_IWOTH);
	string dataset = CONFIG.read<string> ("dataset");
	string genome_filename = input_dir + CONFIG.read<string> ("reference_genome");
	string reads_filename1 = input_dir + CONFIG.read<string> (dataset + "_1");
	string reads_filename2 = input_dir + CONFIG.read<string> (dataset + "_2");
	checkFileExistenceFATAL(genome_filename);
	checkFileExistenceFATAL(reads_filename1);
	checkFileExistenceFATAL(reads_filename2);

	int insert_size = CONFIG.read<int> (dataset + "_IS");
	int dataset_len = CONFIG.read<int> (dataset + "_LEN");
	bool paired_mode = CONFIG.read<bool> ("paired_mode");

	// typedefs :)
	typedef MateReader<Read, ireadstream>::type ReadStream;
	typedef PairedReader<ireadstream> PairedReadStream;
	typedef RCReaderWrapper<PairedReadStream, PairedRead> RCStream;

	// read data ('reads')
	const string reads[2] = { reads_filename1, reads_filename2 };
	ReadStream reader(reads);
	PairedReadStream pairStream(reader, insert_size);
	RCStream rcStream(pairStream);

	// read data ('genome')
	std::string genome;
	{
		ireadstream genome_stream(genome_filename);
		Read full_genome;
		genome_stream >> full_genome;
		genome = full_genome.getSequenceString().substr(0, dataset_len); // cropped
	}
	INFO("Working with " << dataset << " dataset");
	// assemble
	if (paired_mode) {
		debruijn_graph::DeBruijnGraphWithPairedInfoTool<K, RCStream>(rcStream,
				genome, output_dir);
	} else {
		debruijn_graph::DeBruijnGraphTool<K, RCStream>(rcStream, genome,
				output_dir);
	}

	// OK
	return 0;
}


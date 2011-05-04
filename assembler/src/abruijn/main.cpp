#include "toyexamples.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"
#include "libs/getopt_pp/getopt_pp_standalone.h"
#include <iostream>

using namespace GetOpt;
using namespace std;

LOGGER("a");

int main(int argc, char* argv[]) {
	GetOpt_pp options(argc, argv, Include_Environment);

	bool help;
	options >> OptionPresent('h', "help", help);
	if (help) {
		_default_logger->isInfoEnabled();
		_default_logger->setLevel(log4cxx::Level::getError());
		_default_logger->removeAllAppenders();
		_default_logger->closeNestedAppenders();
		_default_logger->warn("\n");
		INFO("\n");
	}
	INFO("Hello, A Bruijn!");
	DEBUG("Using parameters:");

	stringstream usage;
	usage << "\nSet option values like this: \"--key1 value1 --key2 value2 ...\".\n";
	usage << "\nSupported options are:\n";

	int take;
	options >> Option('t', "take", take, 1);
	DEBUG("take = " << take);
	usage << "--take how many minimizers to take from each read (default = 1)\n";

	bool output_single;
	options >> OptionPresent('s', "single", output_single);
	DEBUG("single = " << (output_single ? "true" : "false"));
	usage << "--single (with no value!) whether to visualize vertices non-RC-paired\n";

	int cut;
	options >> Option('c', "cut", cut, -1);
	DEBUG("cut = " << cut);
	usage << "--cut how many first reads from the input file to use (default = -1 = use all reads)\n";

	int mode;
	options >> Option('m', "mode", mode, 0);
	DEBUG("mode = " << mode);
	usage << "--mode bit mask:\n";
	usage << "       1 whether to find minimizers, not local minimizers\n";
	usage << "       2 whether to ensure at least 2 minimizers in each read\n";

	string input_file;
	options >> Option('i', "input", input_file, "");
	string input_file2 = input_file;
	for (int i = input_file2.length() - 1; i >= 0; i--) {
		if (input_file2[i] == '1') {
			input_file2[i] = '2';
			break;
		}
	}
	DEBUG("input = " << input_file << " and " << input_file2);
	usage << "--input first data file. Second file name is obtained by replacing last \"1\" with \"2\"\n";

	string output_file;
	options >> Option('o', "output", output_file, "");
	DEBUG("output = " << output_file);
	usage << "--output output file name";

	if (input_file == "") {
		help = true;
	}
	if (help) {
		std::cout << usage.str();
		return 0;
	}

	std::string file_names[2] = {input_file, input_file2};
	StrobeReader<2, Read, ireadstream> sr(file_names);
	PairedReader<ireadstream> paired_stream(sr, 220);
	SimpleReaderWrapper<PairedReader<ireadstream> > srw(paired_stream);

	abruijn::GraphBuilder gb(srw, take, mode);
	gb.build();

	INFO("Getting statistics...");
	gb.graph.stats();

	INFO("Outputting graph to " << output_file);
	ofstream output_stream(output_file.c_str(), ios::out);
	gb.graph.output(output_stream, !output_single);
	INFO("Done.");
	output_stream.close();

	//ABruijnGraphWithGraphVisualizer ( "ATGTGTGACTTTGTATCGTATTGCGGGCGGCGCGCTTATTGTATGCGTAAATTTGGGTCATATTGATCGTAAAATGCGTATGATGCACTGCA", 6, 3 );
	return 0;
}

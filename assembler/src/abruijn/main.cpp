#include "toyexamples.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"
#include "omni_tools.hpp"
#include "omnigraph.hpp"
#include "visualization_utils.hpp"
#include "libs/getopt_pp/getopt_pp_standalone.h"
#include <iostream>

using namespace GetOpt;
using namespace std;

LOGGER("a");

int main(int argc, char* argv[]) {
	GetOpt_pp options(argc, argv, Include_Environment);
	bool help;
	options >> OptionPresent('h', "help", help);
	stringstream usage;
	usage << "\nSet option values like this: \"--key1 value1 --key2 value2 ...\".\n";
	usage << "\nSupported options are:\n";

	string input_file;
	options >> Option('i', "input", input_file, "");
	string input_file2 = input_file;
	for (int i = input_file2.length() - 1; i >= 0; i--) {
		if (input_file2[i] == '1') {
			input_file2[i] = '2';
			break;
		}
	}
	if (input_file == "") {
		help = true;
	}
	usage << "--input first data file. Second file name is obtained by replacing last \"1\" with \"2\"\n";
	if (!help) {
		INFO("Hello, A Bruijn!");
		INFO("Using parameters:");
		INFO("input = " << input_file << " and " << input_file2);
	}

	string output_file;
	options >> Option('o', "output", output_file, "");
	usage << "--output output file name\n";
	if (!help) INFO("output = " << output_file);

	int take;
	options >> Option('t', "take", take, 1);
	usage << "--take how many minimizers to take from each read (default = 1)\n";
	if (!help) INFO("take = " << take);

	bool output_single;
	options >> OptionPresent('s', "single", output_single);
	usage << "--single (with no value!) whether to visualize vertices non-RC-paired\n";
	if (!help) INFO("single = " << (output_single ? "true" : "false"));

	int cut;
	options >> Option('c', "cut", cut, -1);
	usage << "--cut how many first reads from the input file to use (default = -1 = use all reads)\n";
	if (!help) INFO("cut = " << cut);

	int mode;
	options >> Option('m', "mode", mode, 0);
	usage << "--mode bit mask (default = 0):\n";
	usage << "       1 whether to find minimizers, not local minimizers\n";
	usage << "       2 whether to ensure at least 2 minimizers in each read\n";
	if (!help) {
		if (mode & 1) {
			INFO("mode = find minimizers");
		}
		else {
			INFO("mode = find local minimizers");
		}
		if (mode & 2)
			INFO("mode = ensure at least two minimizers in each read");
	}

	if (help) {
		std::cout << usage.str();
		return 0;
	}

	std::string file_names[2] = {input_file, input_file2};
	StrobeReader<2, Read, ireadstream> sr(file_names);
	PairedReader<ireadstream> paired_stream(sr, 220);
	SimpleReaderWrapper<PairedReader<ireadstream> > srw(paired_stream);
	CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream> > > cr(srw, cut);

	abruijn::GraphBuilderMaster<CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream>>>> gbm(cr, take, mode);
	gbm.build();

//	INFO("Spelling the reference genome");
//	gbm.SpellGenomeThroughGraph(cut + 219);

	INFO("===== Compressing... =====");
	omnigraph::Compressor<omnigraph::Omnigraph> compressor(*gbm.graph());
	compressor.CompressAllVertices();
	INFO(gbm.graph()->size() << " vertices");

//	INFO("===== Getting statistics... =====");
//	gbm.graph()->stats(); TODO

	INFO("Outputting graph to " << output_file);
	ofstream output_stream(output_file.c_str(), ios::out);
//	gbm.graph()->output(output_stream, !output_single); TODO
	gvis::DotPairedGraphPrinter<omnigraph::Omnigraph> printer(*gbm.graph(), "earmarked", output_stream);
	omnigraph::SimpleGraphVisualizer<omnigraph::Omnigraph> sgv(*gbm.graph(), printer);
	sgv.Visualize();
	INFO("Done.");
	output_stream.close();

	//ABruijnGraphWithGraphVisualizer ( "ATGTGTGACTTTGTATCGTATTGCGGGCGGCGCGCTTATTGTATGCGTAAATTTGGGTCATATTGATCGTAAAATGCGTATGATGCACTGCA", 6, 3 );
	return 0;
}

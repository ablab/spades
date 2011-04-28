#include "toyexamples.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"
#include "libs/getopt_pp/getopt_pp_standalone.h"
#include <iostream>

using namespace GetOpt;
using namespace std;

LOGGER("a");

int main(int argc, char* argv[]) {
	INFO("Hello, A Bruijn! Using parameters:");
	GetOpt_pp options(argc, argv, Include_Environment);

	int htake;
	options >> Option('h', "htake", htake, 1);
	INFO("htake = " << htake);

	bool output_single;
	options >> OptionPresent('s', "single", output_single);
	INFO("single = " << (output_single ? "true" : "false"));

	int cut;
	options >> Option('c', "cut", cut, -1);
	INFO("cut = " << cut);

	int mode;
	options >> Option('m', "mode", mode, 0);
	INFO("mode = " << mode);

	string input_file;
	options >> Option('i', "input", input_file);
	string input_file2 = input_file;
	for (int i = input_file2.length() - 1; i >= 0; i--) {
		if (input_file2[i] == '1') {
			input_file2[i] = '2';
			break;
		}
	}
	INFO("input = " << input_file << " and " << input_file2);

	string output_file;
	options >> Option('o', "output", output_file);
	INFO("output = " << output_file);

	std::string file_names[2] = {input_file, input_file2};
	StrobeReader<2, Read, ireadstream> sr(file_names);
	PairedReader<ireadstream> paired_stream(sr, 220);
	SimpleReaderWrapper<PairedReader<ireadstream> > srw(paired_stream);

	abruijn::Graph graph = abruijn::GraphBuilder(srw, htake, mode).build();

	INFO("Getting statistics...");
	graph.stats();

	INFO("Outputting graph to " << output_file);
	ofstream output_stream(output_file.c_str(), ios::out);
	graph.output(output_stream, !output_single);
	output_stream.close();

	//ABruijnGraphWithGraphVisualizer ( "ATGTGTGACTTTGTATCGTATTGCGGGCGGCGCGCTTATTGTATGCGTAAATTTGGGTCATATTGATCGTAAAATGCGTATGATGCACTGCA", 6, 3 );
	return 0;
}

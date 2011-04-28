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

	int label;
	options >> Option('l', "label", label, -1);
	INFO("label = " << label);

	abruijn::Graph graph = abruijn::GraphBuilder(htake).build();

	INFO("Getting statistics...");
	graph.stats();

	INFO("Outputting graph to " << OUTPUT_FILE);
	ofstream outputStream((OUTPUT_FILES + ".dot").c_str(), ios::out);
	graph.output(outputStream);
	outputStream.close();

	//ABruijnGraphWithGraphVisualizer ( "ATGTGTGACTTTGTATCGTATTGCGGGCGGCGCGCTTATTGTATGCGTAAATTTGGGTCATATTGATCGTAAAATGCGTATGATGCACTGCA", 6, 3 );
	return 0;
}

/*
Command line options:
http://code.google.com/p/getoptpp/
http://www.gnu.org/software/hello/manual/libc/Getopt.html
http://sourceforge.net/projects/clp-parser/
*/

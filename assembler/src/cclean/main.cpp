#include <iostream>
#include <map>
#include <exception>

#include "mph_index/kmer_index.hpp"
#include "sequence/seq.hpp"
#include "logger/log_writers.hpp"
#include "memory_limit.hpp"
#include "QcException.hpp"
#include "running_modes.hpp"
#include "ssw/ssw_cpp.h"
#include "config_struct_cclean.hpp"
#include "simple_tools.hpp"

void usage() {
	std::cout << "This tool searches contaminations from UniVec db in provided file with reads" << std::endl;
	std::cout << "Usage: QC-pileline config_path mode:{exact, align, both} UniVec_path Fasta/Fastq.gz" << std::endl;
	std::cout << "Currently only .gz files can be read" << std::endl;
}

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char *argv[]) {

	clock_t start = clock();
	if(5 != argc || (strcmp(argv[2], "exact") && strcmp(argv[2], "align") && strcmp(argv[2], "both"))) {
		usage();
		return 0;
	}

	create_console_logger();
	std::string cfg_filename = argv[1];
	CheckFileExistenceFATAL(cfg_filename);

    INFO("Loading config from " << cfg_filename);
    cclean_cfg::create_instance(cfg_filename);

	const std::string mode(argv[2]);
	std::string db(argv[3]);
	const std::string dt(argv[4]);

	Database * data;
	ireadstream * input;
	try {
		INFO("Reading UniVec db at " << db <<  " ... ");
		data = new Database(db);
		INFO("Done");
		INFO("Init file with reads-to-clean at " << dt << " ... ");
		input = new ireadstream(dt);
		INFO("Done");
	} catch (std::exception& e) {
		ERROR(e.what() << "Make sure that you provided correct path to fasta/fastq file\n");
		return 0;
	}

	INFO("Start matching reads against UniVec ...");
	std::ofstream output(cclean_cfg::get().output_file);
	std::ofstream bed(cclean_cfg::get().output_bed);
	if (!output.is_open() || !bed.is_open()) {
		ERROR("Cannot open output file: " << cclean_cfg::get().output_file << " or " << cclean_cfg::get().output_bed);
		return 0;
	}

	if (!mode.compare("exact")) {
		exactMatch(output, bed, input, data);
	} else if (!mode.compare("align")) {
		alignment(output, bed, input, data);
	} else if (!mode.compare("both")) {
		exactAndAlign(output, bed, input, data);
	}
	output.close();
	bed.close();

	delete data; //NB Don't delete earlier than AhoCorasick, because otherwise all patterns will be deleted
	delete input;

	clock_t ends = clock();
	INFO("Processor Time Spent: " << (double) (ends - start) / CLOCKS_PER_SEC << " seconds.");
	INFO("Goodbye!");
	return 0;
}

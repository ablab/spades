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

//#define TEST

void usage() {
	std::cout << "This tool searches contaminations from UniVec db in provided file with reads" << std::endl;
	std::cout << "Usage: QC-pileline mode:{exact, align, both} UniVec_path Fasta/Fastq.gz" << std::endl;
	std::cout << "Currently only .gz files can be read" << std::endl;
}

void create_console_logger() {
  using namespace logging;

  logger *lg = create_logger("");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

void initSpades() {
	const size_t GB = 1 << 30;
	limit_memory(cfg::get().general_hard_memory_limit * GB);

	// determine quality offset if not specified
	if (!cfg::get().input_qvoffset_opt) {
	  INFO("Trying to determine PHRED offset");
	  int determined_offset = determine_offset(cfg::get().input_single);
	  if (determined_offset < 0) {
		throw QcException("Failed to determine offset! Specify it manually and restart, please!");
	  } else {
		INFO("Determined value is " << determined_offset);
		cfg::get_writable().input_qvoffset = determined_offset;
	  }
	} else {
	  cfg::get_writable().input_qvoffset = *cfg::get().input_qvoffset_opt;
	}
}

void testKmer() {
	try {
    /*    
		KMerData * kmer_data = initSpades();
		KMerDataCounter counter(cfg::get().count_numfiles);
		counter.FillKMerData(*kmer_data);
		std::cout << "total kmers: " << kmer_data->size() << std::endl;
		std::ofstream output("/Users/Kos/Downloads/bio/tmp/kmers.txt");
		kmer_data->binary_write(output);
		output.close();
		std::cout << "Saved" << std::endl;
		delete kmer_data;
    */
	} catch (std::exception& e) {
		ERROR("Caught: "<< e.what());
	}
}

void restoreFromCigar(const std::string& ref, const std::string& query, std::string& out_ref, std::string& out_query, const StripedSmithWaterman::Alignment& a);
void testSSW() {
	  const std::string ref   = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
	  const std::string query = "GCTCTTCCGATCGTGACTGGAGTTCAGACGTGT";

	  // Declares a default Aligner
	  StripedSmithWaterman::Aligner aligner;
	  // Declares a default filter
	  StripedSmithWaterman::Filter filter;
	  // Declares an alignment that stores the result
	  StripedSmithWaterman::Alignment alignment;
	  // Aligns the query to the ref
	  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);
	  std::cout << "===== SSW result =====" << std::endl;
	    std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
	         << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
	         << "Reference start:\t" << alignment.ref_begin << std::endl
	         << "Reference end:\t" << alignment.ref_end << std::endl
	         << "Query start:\t" << alignment.query_begin << std::endl
	         << "Query end:\t" << alignment.query_end << std::endl
	         << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
	         << "Number of mismatches:\t" << alignment.mismatches << std::endl
	         << "Cigar: " << alignment.cigar_string << std::endl;
	    std::cout << "======================" << std::endl;

	    std::cout << ref << std::endl << query << std::endl << std::endl;
	    std::string out_ref, out_query;
	    restoreFromCigar(ref, query, out_ref, out_query, alignment);
	    std::cout << out_ref << std::endl << out_query << std::endl;
}

int main(int argc, char *argv[]) {

	clock_t start = clock();
	create_console_logger();
    INFO("Loading config from " << CONFIG_FILENAME);
    cfg::create_instance(CONFIG_FILENAME);

#ifdef TEST
	testKmer();
//	testSSW();
	return 0;
#endif

	if(4 != argc || (strcmp(argv[1], "exact") && strcmp(argv[1], "align") && strcmp(argv[1], "both"))) {
		usage();
		return 0;
	}

	const std::string mode(argv[1]);
	const std::string dt(argv[3]);
	std::string db(argv[2]);

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
		ERROR(e.what() << "Make sure that you provided correct path to fasta/fastq file in .gz (!)\n");
		return 0;
	}

	INFO("Start matching reads against UniVec ...");
	std::ofstream output(cfg::get().output_file);
	std::ofstream bed(cfg::get().output_bed);
	if (!output.is_open() || !bed.is_open()) {
		ERROR("Cannot open output file: " << cfg::get().output_file << " or " << cfg::get().output_bed);
		return 0;
	}

	if (!strcmp(argv[1], "exact")) {
		exactMatch(output, bed, input, data);
	} else if (!strcmp(argv[1], "align")) {
		alignment(output, bed, input, data);
	} else if (!strcmp(argv[1], "both")) {
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

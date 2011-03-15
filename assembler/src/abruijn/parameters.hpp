#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <sstream>
#include "ireadstream.hpp"

// K-mer size
#define K 25

// ===== Input data =====
#define FILENAMES make_pair("./data/MG1655-K12_emul1.fasta.gz", "./data/MG1655-K12_emul2.fasta.gz")
//#define FILENAMES make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz")
// How many reads we process
#define CUT 200000
// Mate-pair read size
#define MPSIZE 100

// ===== Output =====
#define OUTPUT_FILE std::string("./output/") + inputFileNames.first[7] + itoa(K) + "_" + itoa(CUT)

// ===== Visualization =====
// Define OUTPUT_PAIRED if you want GraphViz to glue together complementary vertices
//#define OUTPUT_PAIRED
// How many nucleotides are put in the vertex label in GraphViz
#define LABEL 5

// ===== Hashing parameters =====
#define HASH_XOR 1845724623
#define HASH_X(v) ((v << 5) - v)
typedef unsigned int hash_t;

// ===== Debug parameters =====
#define VERBOSE(n, message) if ((n & 8191) == 8191) INFO(n << message)

static string itoa(int n) {
	std::stringstream ss;
	ss << n;
	return ss.str();
}

static pair<std::string, std::string> inputFileNames = FILENAMES;
static ireadstream<MPSIZE, 2> irs(inputFileNames.first.c_str(), inputFileNames.second.c_str());

#endif /* PARAMETERS_H_ */

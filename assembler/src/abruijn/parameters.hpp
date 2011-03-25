#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <fstream>
#include <sstream>
#include "ireadstream.hpp"

// K-mer size
#define K 25

// ===== Input data =====
#define FILENAMES make_pair("./data/MG1655-K12_emul1.fastq.gz", "./data/MG1655-K12_emul2.fastq.gz")
//#define FILENAMES make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz")
// How many reads we process
#define CUT 1000000
// Mate-pair read size
#define MPSIZE 100

// ===== Visualization =====
// Define OUTPUT_PAIRED if you want GraphViz to glue together complementary vertices
#define OUTPUT_PAIRED
// How many nucleotides are put in the vertex label in GraphViz
#define LABEL 5

// ===== Output =====
#define OUTPUT_FILE std::string("./output/") + OUTPUT_FILE_NAME
#define OUTPUT_FILE_NAME inputFileNames.first[7] + toString(K) + "_" + toString(CUT) + OUTPUT_FILE_SUFFIX
#ifdef OUTPUT_PAIRED
	#define OUTPUT_FILE_SUFFIX ""
#endif
#ifndef OUTPUT_PAIRED
	#define OUTPUT_FILE_SUFFIX "_s"
#endif

// ===== Hashing parameters =====
#define HASH_XOR 1845724623
#define HASH_X(v) ((v << 5) - v)
typedef unsigned int hash_t;

// ===== Debug parameters =====
#define VERBOSE(n, message) if ((n & 8191) == 8191) INFO(n << message)

template <typename T>
string toString(T t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

static pair<std::string, std::string> inputFileNames = FILENAMES;
static ireadstream<MPSIZE, 2> inputStream(inputFileNames.first.c_str(), inputFileNames.second.c_str());
static ofstream outputStream((OUTPUT_FILE + ".dot").c_str(), ios::out);


#endif /* PARAMETERS_H_ */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <fstream>
#include <sstream>
#include "ireadstream.hpp"

/**
 * K-mer size
 */
#define K 25

// ===== Directories ===== //
#define INPUT_DIRECTORY std::string("./data/input/")
#define OUTPUT_DIRECTORY std::string("./data/abruijn/")

// ===== Input data ===== //
#define INPUT_DATA_SET std::string("MG1655-K12_emul")
//#define INPUT_DATA_SET std::string("s_6_")
/**
 * Mate-pair read size
 */
#define MPSIZE 100
/**
 * How many reads we process
 */
#define CUT 400

// ===== Visualization ===== //
/**
 * Define OUTPUT_PAIRED if you want GraphViz to glue together complementary vertices
 */
#define OUTPUT_PAIRED
/**
 * How many nucleotides are put in the vertex label in GraphViz
 */
//#define LABEL 5

template <typename T>
string toString(T t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

// ===== Naming conventions ===== //
#define INPUT_FILES INPUT_DIRECTORY + INPUT_DATA_SET + "1.fastq.gz", INPUT_DIRECTORY + INPUT_DATA_SET + "2.fastq.gz"
#define OUTPUT_FILE INPUT_DATA_SET + toString(K) + "_" + toString(CUT) + OUTPUT_FILE_SUFFIX
#ifdef OUTPUT_PAIRED
	#define OUTPUT_FILE_SUFFIX ""
#endif
#ifndef OUTPUT_PAIRED
	#define OUTPUT_FILE_SUFFIX "_s"
#endif
#define OUTPUT_FILES OUTPUT_DIRECTORY + OUTPUT_FILE

// ===== Hashing parameters ===== //
#define HASH_XOR 1845724623
#define HASH_X(v) ((v << 5) - v) // = v * 31
// Type of hash values.
typedef unsigned int hash_t;
// Maximum value of type hash_t
const hash_t maxHash = -1;

// ===== Debug parameters ===== //
#define VERBOSE(n, message) if ((n & 8191) == 8191) INFO(n << message)

#endif /* PARAMETERS_H_ */

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
#define INPUT_DIRECTORY std::string("./data/input/cropped/")
#define OUTPUT_DIRECTORY std::string("./data/abruijn/")

// ===== Input data ===== //
//#define INPUT_DATA_SET std::string("MG1655-K12_emul")
//#define INPUT_DATA_SET std::string("s_6_")
#define INPUT_DATA_SET std::string("s_6.first10000_")

// ===== Naming conventions ===== //
#define INPUT_FILES INPUT_DIRECTORY + INPUT_DATA_SET + "1.fastq.gz", INPUT_DIRECTORY + INPUT_DATA_SET + "2.fastq.gz"
#define OUTPUT_FILE INPUT_DATA_SET + toString(K)
#define OUTPUT_FILES OUTPUT_DIRECTORY + OUTPUT_FILE

// ===== Debug parameters ===== //
#define VERBOSE(n, message) if (n % 10000 == 0 && n > 0) INFO(n << message)

#endif /* PARAMETERS_H_ */

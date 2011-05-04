/*
 * config.hpp
 *
 *  Created on: 03.05.2011
 *      Author: vyahhi
 */

#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include <tr1/tuple>
#include <string>

const std::string INPUT_DIR = "./data/input/";
const std::string ECOLI_FILE = "./data/input/MG1655-K12.fasta.gz";
//const std::string READS_FILE = "./data/input/s_6.first100000_1.fastq.gz";
const std::tr1::tuple<std::string, std::string, size_t, int> QUAKE_CROPPED_10_3   = std::tr1::make_tuple<std::string, std::string, size_t, int>(INPUT_DIR + "s_6.first1000_1.fastq.gz",   INPUT_DIR + "s_6.first1000_2.fastq.gz", 220, 1000);
const std::tr1::tuple<std::string, std::string, size_t, int> QUAKE_CROPPED_10_4   = std::tr1::make_tuple<std::string, std::string, size_t, int>(INPUT_DIR + "s_6.first10000_1.fastq.gz",  INPUT_DIR + "s_6.first10000_2.fastq.gz", 220, 10000);
const std::tr1::tuple<std::string, std::string, size_t, int> QUAKE_CROPPED_10_5   = std::tr1::make_tuple<std::string, std::string, size_t, int>(INPUT_DIR + "s_6.first100000_1.fastq.gz", INPUT_DIR + "s_6.first100000_2.fastq.gz", 220, 100000);
const std::tr1::tuple<std::string, std::string, size_t, int> QUAKE_CROPPED_4_10_5 = std::tr1::make_tuple<std::string, std::string, size_t, int>(INPUT_DIR + "s_6.first400000_1.fastq.gz", INPUT_DIR + "s_6.first400000_2.fastq.gz", 220, 4000000);

#define K 27
#define R 100
#define I 220
#define DE_BRUIJN_DATA_FOLDER "./data/debruijn/"
#define INPUT QUAKE_CROPPED_10_5

#endif /* CONFIG_HPP_ */

/*
 * parameters.h
 *
 *  Created on: 22.02.2011
 *      Author: misha
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "seq.hpp"

#define MPSIZE 100
#define K 25
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;
typedef Seq<MPSIZE> Read;

using namespace std;

static pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");

#endif /* PARAMETERS_H_ */

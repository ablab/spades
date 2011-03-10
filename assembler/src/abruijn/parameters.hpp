#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "seq.hpp"
#include "ireadstream.hpp"

#define MPSIZE 100
#define K 25
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;
typedef Seq<MPSIZE> Read;

using namespace std;

static pair<string,string> filenames = make_pair("./data/MG1655-K12_emul1.fasta.gz", "./data/MG1655-K12_emul2.fasta.gz");
//static pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");
//static ireadstream<MPSIZE, 2> irs(filenames.first.c_str(), filenames.second.c_str());

static ireadstream<MPSIZE, 2> datastream() {
	return ireadstream<MPSIZE, 2>(filenames.first.c_str(), filenames.second.c_str());
}

#endif /* PARAMETERS_H_ */

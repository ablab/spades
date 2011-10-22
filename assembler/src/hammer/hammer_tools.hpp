/*
 * hammer_tools.hpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "union.hpp"
#include "sequence/seq.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"
#include "position_kmer.hpp"
#include "subkmers.hpp"

using namespace std;

#define MAX_INT_64 1000000000000000000

#define TIMEDLN(a) print_time(); cout << a << endl

double oct2phred(string qoct, int qvoffset);
string encode3toabyte (const string & s);
void print_time();

/// join two maps
void join_maps(KMerStatMap & v1, const KMerStatMap & v2);

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const PositionRead &r, hint_t readno, KMerStatMap *v);

void AddKMerNos(const PositionRead &r, hint_t readno, vector<KMerNo> *v);

void DoPreprocessing(int tau, string readsFilename, int nthreads, vector<KMerCount*> * kmers, KMerNoHashMap * km);
void DoSplitAndSort(int tau, int nthreads, vector< vector<hint_t> > * vs, vector<KMerCount> * kmers, vector<SubKMerPQ> * vskpq);
void outputReads(bool paired, const char * fname, const char * fname_bad, const char * fname_right = NULL, const char * fname_right_bad = NULL, const char * fname_left_unpaired = NULL, const char * fname_right_unpaired = NULL);

/**
  * correct a read in place
  * @return how many nucleotides have been changed
  */
size_t CorrectRead(const KMerNoHashMap & km, const vector<KMerCount*> & kmers, hint_t readno, Read & r, bool & isGood, ofstream * ofs = NULL);
bool internalCorrectReadProcedure( const Read & r, const hint_t readno, const string & seq,
		const vector<KMerCount*> & km, const PositionKMer & kmer, const uint32_t pos, const KMerStat & stat,
		vector< vector<int> > & v, int & left, int & right, bool & isGood, ofstream * ofs );

/**
  * make a step of iterative reconstruction
  * @return number of new solid k-mers
  */
size_t IterativeReconstructionStep(int nthreads, const vector<KMerCount*> & kmers, ostream * ofs = NULL);

/**
 * This function reads reads from the stream and splits them into
 * k-mers. Then k-mers are written to several file almost
 * uniformly. It is guaranteed that the same k-mers are written to the
 * same files.
 * Different from quake_count in that it works with position reads
 * @param dirprefix where to put the temporary files
 * @param iter_count no. of current iteration
 */
void SplitToFiles(string dirprefix, int iter_count);

/**
 * process a single file with kmers divided by hashes
 * output results into kmerno_file
 */
void ProcessKmerHashFile( ifstream * inStream, ofstream * kmerno_file, 	hint_t & kmer_num );

/**
 * fill in kmerno vector
 */
void fillInKmersFromFile( const string & fname, vector<hint_t> *kmernos );

string getFilename( const string & dirprefix, const string & suffix );
string getFilename( const string & dirprefix, int iter_count, const string & suffix );
string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num );

void getGlobalConfigParameters( const string & config_file );

void dumpBlob( const string & fname );
void loadBlob( const string & fname );

#endif


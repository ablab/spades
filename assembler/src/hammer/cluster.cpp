/*
 * cluster.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#include<cassert>
#include<cmath>
#include<string.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<fcntl.h>

using namespace std;
ostringstream oMsg;
string sbuf;
#include "defs.hpp"
#include "union.hpp"
#include "hammer_config.hpp"

long get_filesize(string filename) {
	int file=0;
	if((file=open(filename.c_str(),O_RDONLY)) < -1) {
		cerr << "File not found: " << filename;
		exit(1);
	}
	struct stat fileStat;
	if(fstat(file,&fileStat) < 0)  {
		cerr << "Can't fstat file: " << filename;
		exit(1);
	}
	return fileStat.st_size;
}

long load_file_to_memory(const char *filename, char **result) 
{ 
	long size = get_filesize(filename);
	FILE *f = fopen(filename, "rb");
	if (f == NULL) { 
		*result = NULL;
		cout << "Failed opening file: " << filename << endl;
		exit(1);
	} 
	*result = (char *)malloc(size+1);
	if (size != fread(*result, sizeof(char), size, f)) { 
		cout << "Failed reading file: " << filename << endl;
		free(*result);
		exit(1);
	} 
	fclose(f);
	(*result)[size] = 0;
	return size;
}

int tau;
int nthreads;
long numDistinctKmers;
string kmerFile;
string splitsBase;
int kmerfileLinelen;
char * reads;

class entryClass { 
public:
	string key;
	int id;
	string seq;
	string count;
	string freq;
};


void get_full_kmer(entryClass & x) {
	char s[K+1];
	long offset = x.id * long(kmerfileLinelen);
	strncpy(s, &reads[offset], K);
	s[K] = 0;
	x.seq = s;
	strncpy(s, &reads[offset + K], 6);
	s[6] = 0;
	x.count = s;
	strncpy(s, &reads[offset + K + 6], 9);
	s[9] = 0;
	x.freq = s;
	// cout << offset << "  " << x.seq << ":" << x.count << ":" << x.freq << endl;
	return;
}

void process_block(unionFindClass * uf, vector<entryClass> & block) {
	//cout << "processing block of size " << block.size() << "\n";
	for (int i = 0; i < block.size(); i++) {
		get_full_kmer(block[i]);
	}
	for (int i = 0; i < block.size(); i++) {
		uf->find_set(block[i].id);
		for (int j = i + 1; j < block.size(); j++) {
			//cout << "comparing " << block[i].seq << " and " << block[j].seq << "\n";
			if (hamdist(block[i].seq, block[j].seq, SAME_STRAND, tau ) <= tau) {
				uf->unionn(block[i].id, block[j].id);
			}
		}
	}
	return;
}
typedef pair<int, unionFindClass *> paramType;

void * onethread(void * params) {
	int thread = ((paramType *) params)->first;
	unionFindClass * uf = new unionFindClass(numDistinctKmers);
	((paramType *) params)->second = uf;
	char threadLabel[5];
	sprintf(threadLabel, "%02d", thread);

	cerr << "Processing split files (" << thread << ")...\n";
	ifstream inf;
	// cout << splitsBase + "." + threadLabel << "\n";
	open_file(inf, splitsBase + "." + threadLabel);

	string sbuf;
	entryClass last;
	vector<entryClass> block;
	int counter=0;
	while (getline(inf, sbuf)) { 
		if (++counter % 1000000 == 0) cerr << "Processed (" << thread << ") " << add_commas(counter) << ", ";
		istringstream line(sbuf);
		entryClass cur;
		line >> cur.key >> cur.id;
		//cout << cur.key << " " << cur.id << " " << cur.id * long(kmerfileLinelen) << "\n";
		if (last.key == cur.key) { //add to current reads
			block.push_back(cur);
		} else {
			process_block(uf, block);
			block.clear();
			block.push_back(cur);
		}
		last = cur;
	}
	process_block(uf, block);
	inf.close();
	cerr << "Finished(" << thread << ") " << endl;
	pthread_exit(NULL);
}

void merge(vector<paramType> params, string prefix) {
	cerr << "Merging union find files...\n";
	unionFindClass * ufMaster = new unionFindClass(numDistinctKmers);
	vector<string> row;
	vector<vector<int> > classes;
	for (int i = 0; i < nthreads; i++) {
		classes.clear();
		unionFindClass * ufSlave = params[i].second;
		ufSlave->get_classes(classes);
		delete ufSlave;
		for (int j = 0; j < classes.size(); j++) {
			int first = classes[j][0];
			ufMaster->find_set(first);
			for (int k = 0; k < classes[j].size(); k++) {
				ufMaster->unionn(first, classes[j][k]);
			}
		}
	}

	cerr << "Outputting clusters...\n";
	classes.clear();
	ufMaster->get_classes(classes);
	delete ufMaster;
	entryClass x;
	ofstream outf;
	open_file(outf, prefix + "reads.uf");
	for (int i = 0; i < classes.size(); i++) {
		for (int j = 0; j < classes[i].size(); j++) {
			x.id = classes[i][j];
			get_full_kmer(x);
			outf << "ITEM\t" << i << "\t" << x.id << "\t" << x.seq << "\t" << x.count << "\t" << x.freq <<  endl;
		}
		outf << endl;
	}
	outf.close();
}

int main(int argc, char * argv[]) {
	tau = atoi(argv[1]);
	kmerFile = argv[2];
	splitsBase = argv[3];
	long memlim = atol(argv[4]);
	string mergeOnly = argv[5];
	kmerfileLinelen = K + 16;
	nthreads = atoi(argv[6]);
	string dirprefix = argv[7];
	//int nthreads = tau + 1;

	//this depends on the implementation of UnionFindClass
	long l = get_filesize(kmerFile);
	cout << l << " " << kmerfileLinelen << "\n";
	assert (l % kmerfileLinelen == 0);
	numDistinctKmers = l / kmerfileLinelen;
	cout << numDistinctKmers << "\n";
	//cerr << "numDistinct = " << numDistinctKmers << endl;

	//load reads
	cerr << "Reading in kmers...\n";
	load_file_to_memory(kmerFile.c_str(), &reads);
	// cout << "\n\n" << reads << "\n\n";

	if (mergeOnly == "1") {
		cerr << argv[0] << ": mergeOnly mode not supported, exiting.\n";
		exit(1);
	}

	long maxKmers =  (memlim - l - 2*sizeof(int)) / (nthreads * 2*sizeof(int));  //accounting for union find 
	cerr << argv[0] << ": maxKmers = " << add_commas(maxKmers) << ", numDistinctKmers = " << add_commas(numDistinctKmers) << endl;
	if (numDistinctKmers > maxKmers) {
		cerr << "Won't fit into memory, don't even try fool.\n";
		exit(1);
	}

	//start threads
	vector<pthread_t> thread(nthreads);
	vector<paramType> params(nthreads);
	for (int i = 0; i < nthreads; i++) {
		cout << i << " out of " << nthreads << "\n";
		params[i].first = i;
		pthread_create(&thread[i], NULL, onethread, (void *) &params[i]);
	}
	cout << "Initiated all threads.\n";
	for (int i = 0; i < nthreads; i++) {
		pthread_join(thread[i], NULL);
	}

	cout << "printing to " << dirprefix.data() << "/reads.uf\n";
	merge(params, dirprefix + "/");

	delete reads;
	return 0;
}



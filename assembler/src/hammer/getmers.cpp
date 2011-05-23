/*
 * getmers.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
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

#include "defs.hpp"
#include "union.hpp"
#include "hammer_config.hpp"

using namespace std;
using namespace __gnu_cxx;
ostringstream oMsg;
string sbuf;

int kmersize = 0;
char * reads;

struct Kmer {
	long loc;
	int mult;
	float freq;
};


int comp (const void * k1, const void * k2) {
	long x = ((Kmer *) k1)->loc;
	long y = ((Kmer *) k2)->loc;

	if (x==y) return 0;
	bool flipx = false;  // are we revcomping x?
	bool flipy = false;
	bool detx = false; //have we determine whether we are revcomping x?
	bool dety = false;
	for (int i = 0; i < kmersize; i++) {
		char xchar, ychar;
		if (detx && !flipx) {
			xchar = reads[x+i];
		} else if (detx && flipx) {
			xchar = revcomp(reads[x + kmersize - i - 1]);
		} else { //!detx
			if (reads[x+i] == revcomp(reads[x+kmersize -i - 1])) {
				xchar = reads[x+i];
			} else {
				detx = true;
				flipx = (reads[x+i] > revcomp(reads[x+kmersize -i - 1]));
				if (flipx) {
					xchar = revcomp(reads[x+kmersize -i - 1]);
				} else {
					xchar = reads[x + i];
				}
			}
		}
		if (dety && !flipy) {
			ychar = reads[y+i];
		} else if (dety && flipy) {
			ychar = revcomp(reads[y + kmersize - i - 1]);
		} else { //!dety
			if (reads[y+i] == revcomp(reads[y+kmersize -i - 1])) {
				ychar = reads[y+i];
			} else {
				dety = true;
				flipy = (reads[y+i] > revcomp(reads[y+kmersize -i - 1]));
				if (flipy) {
					ychar = revcomp(reads[y+kmersize -i - 1]);
				} else {
					ychar = reads[y + i];
				}
			}
		}

		if (xchar < ychar) return -1;
		if (xchar > ychar) return 1;
	}
	return 0;
}



int main(int argc, char * argv[]) {
	kmersize = atoi(argv[1]);
	unsigned long readslim = strtoul(argv[2], NULL, 0);
	unsigned long totlim = strtoul(argv[3], NULL, 0);
	string preproFilename = argv[4];

	cerr << "Reading in prepro file...\n" << sizeof(Kmer) << "\n";
	
	ireadstream ifs(readsFilename.data());
	ofstream ofs;
	ofs.open(outFilename.data());
	Read r;
	while (!ifs.eof()) {
		ifs >> r;
		string corrs = correctSequence(r.getSequence(), ufc_map);
		ofs << "@" << r.getName() << endl << corrs.data() << endl << "+" << endl << r.getPhredQualityString() << endl;
	}

	
	KMerHashMap ufc_map;
	UFStream ufs(ufFilename);
	MyUFC cur_ufc;
	while(!ufs.eof()) {
		ufs >> cur_ufc;
		for (size_t i=0; i<cur_ufc.size(); ++i) {
			if (i != cur_ufc.center())
				ufc_map.insert(make_pair(cur_ufc.seq(i), cur_ufc.seq(cur_ufc.center())));
		}
	}	
	cout << ufc_map.size() << "\n";

	
	long maxKmers = (totlim - readslim) / sizeof(Kmer);
	cerr << "maxKmers = " << totlim << " "<< readslim << " " << (maxKmers) << endl;
	Kmer * indices = new Kmer[maxKmers];
	reads = new char[readslim]; 
	vector<string> row;
	long counter = 0;
	Kmer kmer;
	long counter2 = 0;
	long indicesSize = 0;
	while (get_row_whitespace(cin, row)) {
		counter2++;
		memcpy(&reads[counter], row[0].data(), row[0].length());
		for (int i = 0; i < max(0, (int) row[0].length() - kmersize + 1); i++) { //don't know why the max is necessary, but it wouldn't work otherwise
			kmer.loc = counter + i;
			kmer.mult = 1;
			kmer.freq = atof(row[1].c_str());
			if (indicesSize > maxKmers) {
				cerr << "Out of memory, after " << counter2 << " lines of your reads.prepro.\n";
				exit(1);
			}
			indices[indicesSize] = kmer;
			indicesSize++;
		}
		counter += row[0].length();
	}
	
	cerr << "Sorting kmers...\n";
	qsort(indices, indicesSize, sizeof(Kmer), comp);

	cerr << "Removing duplicates...\n";
	long upos = 0;
	for (long i = 1; i < indicesSize; i++) {
		if (0 == comp(&indices[upos], &indices[i])) {
			indices[upos].mult += indices[i].mult;
			indices[upos].freq += indices[i].freq;
		} else {
			upos++;
			indices[upos] = indices[i];
		}
	}
	indicesSize = upos + 1;

	cerr << "Outputing kmers file...\n";
	for (long i = 0; i < indicesSize; i++) {
		string seq = "";
		for (int j = 0; j < kmersize; j++) {
			seq += reads[indices[i].loc + j];
		}
		assert(indices[i].mult < 100000);
		assert(indices[i].freq < 100000);
		printf("%s %5i %8.2f\n", rcnorm(seq).c_str(), indices[i].mult, indices[i].freq);
	}
	delete indices;
}




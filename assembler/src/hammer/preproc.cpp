/*
 * preproc.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#include<cmath>
#include<string>
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
#include<cassert>

#include "hammer_config.hpp"
#include "defs.hpp"
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "sequence/seq.hpp"

using namespace std;
ostringstream oMsg;

int qvoffset;


double oct2phred(string qoct)  {
	float freq = 1;	
	for (size_t i = 0; i < qoct.length(); i++) {
		freq *= 1 - pow(10, -float(qoct[i] - qvoffset)/10.0);
	}

	return freq;
}               

string encode3toabyte (const string & s)  {
	string retval;
	char c = 48;
	int weight = 16;
	size_t i;
	for (i = 0; i < s.length(); i += 1) {
		if (i % 3 == 0) {
			c= 48;
			weight = 16;
		}
		c += weight * nt2num(s[i]);
		weight /= 4;
		if (i % 3 == 2) retval += c;
	}
	if (i % 3 != 0) retval += c;
	return retval;
}


/**
 * add k-mers from read to map
 */
void addKMers(const Read & r, KMerStatVector & v) {
	size_t i = 0;
	KMerStatMap::iterator it;
	float freq = oct2phred(r.getPhredQualityString(qvoffset));
	string s = r.getSequenceString();
	while (true) {
		i = r.firstValidKmer(i, K);
		if (i+K > r.size()) break;
		KMer kmer = KMer(r.getSubSequence(i, K));
		while (true) {
			/*it = m.find(kmer);
			if (it != m.end()) {
				it->second.count++;
				it->second.freq += freq;
			} else {
				pair<KMer, KMerStat> p;
				p.first = kmer;
				p.second.count = 1; p.second.freq = freq;
				m.insert(p);
			}*/
			KMerStat stat; stat.count = 1; stat.freq = freq;
			v.push_back( make_pair(kmer, stat) );
			if (i+K < r.size() && is_nucl(s[i+K])) {
				kmer = kmer << r[i+K];
				++i;
			} else {
				i = i+K;
				break;
			}
		}
	}
}

int main(int argc, char * argv[]) {

	//cout << Reverse(Complement("AAAAAAAAAAAAAAAATCGGAGGAATAAAAATAGAATAGTATACAAAGGCCCCTT")) << "\n";
	//return 0;
	//int starttrim = atoi(argv[1]); 
	//int endtrim = atoi(argv[2]);  //number of characters to keep
	int tau = atoi(argv[1]);
	qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string outFilename = argv[4];
	string kmerFilename = argv[5];
	
	KMerStatVector v;

	ireadstream ifs(readsFilename.data(), qvoffset);
	ofstream ofs;
	ofs.open(outFilename.data());
	Read r;
	while (!ifs.eof()) {
		ifs >> r;
		// trim the reads for bad quality and process only the ones with at least K "reasonable" elements
		if (r.trimBadQuality() >= K) {
			// cout << r.getSequenceString() << "\n" << (!r).getSequenceString() << "\n\n";
			addKMers(r, v);
			addKMers(!r, v);
			ofs << "@" << r.getName() << endl << r.getSequenceString().data() << endl << "+" << endl << r.getPhredQualityString(qvoffset) << endl;
		}
	}
	ifs.close();
	ofs.close();
	cout << "All k-mers added to the vector.\n";
	
	sort(v.begin(), v.end(), KCgreater);
	size_t j = 0;
	for (size_t i=1; i<v.size(); ++i) {
		if ( v[i].first == v[j].first ) {
			v[j].second.count += v[i].second.count;
			v[j].second.freq += v[i].second.freq;
		} else {
			++j;
		}
	}
	v.resize(j+1);
	cout << "All k-mers sorted.\n";	
	
	FILE* f = fopen(kmerFilename.data(), "w");
	size_t counter = 0;
	vector<StringCountVector> vs(tau+1);
	for (KMerStatVector::const_iterator it = v.begin(); it != v.end(); ++it) {
		fprintf(f, "%s %5u %8.2f\n", it->first.str().data(), it->second.count, it->second.freq);
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += nucl(it->first[i]);
			}
			// cout << it->first.str().data() << ", " << j << ": " << sub << "   " << encode3toabyte(sub) << "\n";
			vs[j].push_back(make_pair(encode3toabyte(sub), counter));
		}
		++counter;
	}
	fclose(f);
	cout << "Auxiliary vectors loaded.\n";	
	
	for (int j=0; j<tau+1; ++j) {
		sort(vs[j].begin(), vs[j].end(), SCgreater);
	}	
	cout << "Auxiliary vectors sorted.\n";	
		
	for (int j=0; j<tau+1; ++j) {
		stringstream fname;
		fname << kmerFilename << "." << j;
		ofs.open(fname.str().data());
		for (StringCountVector::const_iterator it = vs[j].begin(); it != vs[j].end(); ++it) {		
			ofs << it->first << "\t" << it->second << endl;
		}
		ofs.close();
	}
	
	
	
	return 0;
}



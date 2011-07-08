/*
 * preproc.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#include<omp.h>
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
void addKMers(const Read & r, KMerStatMap & v) {
	KMerStatMap::iterator it;
	float freq = oct2phred(r.getPhredQualityString(qvoffset));
	string s = r.getSequenceString();
	size_t i=0;
	while (true) {
		i = r.firstValidKmer(i, K);
		if (i+K > r.size()) break;
		KMer kmer = KMer(r.getSubSequence(i, K));
		while (true) {
			it = v.find(kmer);
			if (it != v.end()) {
				it->second.count++;
				it->second.freq += freq;
			} else {
				pair<KMer, KMerStat> p;
				p.first = kmer;
				p.second.count = 1; p.second.freq = freq;
				v.insert(p);
			}
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

void join_maps(KMerStatMap & v1, const KMerStatMap & v2) {
	KMerStatMap::iterator itf;
	for (KMerStatMap::const_iterator it = v2.begin(); it != v2.end(); ++it) {
		itf = v1.find(it->first);
		if (itf != v1.end()) {
			itf->second.count = itf->second.count + it->second.count;
			itf->second.freq  = itf->second.freq  + it->second.freq;
		} else {
			pair<KMer, KMerStat> p;
			p.first = it->first;
			p.second.count = it->second.count;
			p.second.freq = it->second.freq;
			v1.insert(p);
		}
	}
}

#define READ_BATCH_SIZE 1000000
#define BATCHES_PER_MAP 4
#define MAX_INT_64 15000000000000000000


class ReadStatMapContainer {
public:
	ReadStatMapContainer(const vector<KMerStatMap> & vv) : v_(vv) { init(); }

	void init() {
		i_.clear();
		for (size_t i=0; i<v_.size(); ++i) {
			i_.push_back(v_[i].begin());
		}

	}

	pair<KMer, KMerStat> next() {
		pair<KMer, KMerStat> p;
		KMerStatMap::const_iterator imin = cur_min();
		if (imin == v_[0].end()) {
			p.second.count = MAX_INT_64;
			return p;
		}
		p.first = imin->first; p.second.count = 0; p.second.freq = 0;
		//cout << "Next. Min=" << imin->first.str().data() << "\n";
		for (size_t i=0; i<v_.size(); ++i) {
			/*if (i_[i] != v_[i].end())
				cout << "   " << i << " = " << i_[i]->first.str().data() << "\n";*/
			if (i_[i]->first == p.first) {
				//cout << i << "=" << i_[i]->second.count << " ";
				p.second.count += i_[i]->second.count;
				p.second.freq  += i_[i]->second.freq;
				++i_[i];
			}
		}
		//cout << " = " << p.second.count << "\n\n";
		return p;
	}

private:
	const vector<KMerStatMap> & v_;
	vector<KMerStatMap::const_iterator> i_;

	const KMerStatMap::const_iterator & cur_min() {
		int min=0;
		for (size_t i=1; i<v_.size(); ++i) {
			if (i_[i] != v_[i].end()) {
				if (i_[min] == v_[min].end() || KMerLess(i_[i]->first, i_[min]->first)) {
					min = i;
				}
			}
		}
		return i_[min];
	}

};


int main(int argc, char * argv[]) {

	int tau = atoi(argv[1]);
	qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string outFilename = argv[4];
	string kmerFilename = argv[5];
	int nthreads = atoi(argv[6]);

	/*
	
	vector<KMerStatMap> vv;
	for (int i=0; i<nthreads; ++i) {
		KMerStatMap v; v.clear();
		vv.push_back(v);
	}

	cout << "Starting preproc.\n";
	ireadstream ifs(readsFilename.data(), qvoffset);
	ofstream ofs;
	ofs.open(outFilename.data());
	Read r;
	size_t tmpc = 0;
	size_t cur_maps = 0;
	vector<Read> rv;
	while (!ifs.eof()) {
		// reading a batch of reads
		for (int thr = 0; thr < READ_BATCH_SIZE; ++thr) {
			ifs >> r; 
			if (r.trimBadQuality() >= K) {
				rv.push_back(r);
				ofs << "@" << r.getName() << endl << r.getSequenceString().data() << endl << "+" << endl << r.getPhredQualityString(qvoffset) << endl;
			}
			if (ifs.eof()) break;
		}

		// trim the reads for bad quality and process only the ones with at least K "reasonable" elements
		// we now do it in parallel

		++tmpc;
		cout << "Batch " << tmpc << " read.\n"; flush(cout);
		#pragma omp parallel for shared(rv, vv, ofs) num_threads(nthreads)
		for(int i=0; i<rv.size(); ++i) {
			addKMers(rv[i], vv[omp_get_thread_num() + cur_maps * nthreads]);
			addKMers(!(rv[i]), vv[omp_get_thread_num() + cur_maps * nthreads]);
		}
		cout << "Batch " << tmpc << " added.\n"; flush(cout);
		rv.clear();
	}
	ifs.close();
	ofs.close();
	cout << "All k-mers added to maps.\n"; flush(cout);

	ReadStatMapContainer rsmc(vv);
	for (int i=0; i<vv.size(); ++i) cout << "size(" << i << ")=" << vv[i].size() << "\n"; flush(cout);

	FILE* f = fopen(kmerFilename.data(), "w");
	size_t counter = 0;
	vector<StringCountVector> vs(tau+1);
	//for (KMerStatMap::const_iterator it = vv[0].begin(); it != vv[0].end(); ++it) {
	for (pair<KMer, KMerStat> p = rsmc.next(); p.second.count < MAX_INT_64; p = rsmc.next()) {
		fprintf(f, "%s %5u %8.2f\n", p.first.str().data(), p.second.count, p.second.freq);
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += nucl(p.first[i]);
			}
			vs[j].push_back(make_pair(encode3toabyte(sub), counter));
		}
		++counter;
	}
	fclose(f);  */
	

	int effective_threads = min(nthreads, tau+1);	

	FILE* f = fopen(kmerFilename.data(), "r");
	size_t counter = 0;
	vector<StringCountVector> vs(tau+1);
	char kmerstr[K], myline[K + 30];
	int count; float freq;
	while (!feof(f)) {
		fgets(myline, K+30, f); if (feof(f)) break;
		sscanf(myline, "%s %5u %8.2f", kmerstr, &count, &freq);
		
		#pragma omp parallel for shared(vs) num_threads(effective_threads)
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += kmerstr[i];
			}
			vs[j].push_back(make_pair(encode3toabyte(sub), counter));
		}
		++counter;
		if (counter % 1000000 == 0) cout << "Read " << counter << " kmers.\n"; flush(cout); 
	}
	fclose(f);
	cout << "Auxiliary vectors loaded from file " << kmerFilename.data() << ".\n"; flush(cout);

		
	#pragma omp parallel for shared(vs) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		sort(vs[j].begin(), vs[j].end(), SCgreater);
	}
	
	cout << "Auxiliary vectors sorted.\n"; flush(cout);
	
	ofstream ofs;

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



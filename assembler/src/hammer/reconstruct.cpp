/*
 * center.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#include<cassert>
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
#include<stdexcept>

#include "hammer_config.hpp"
#include "defs.hpp"
#include "hammerread.hpp"
#include "mathfunctions.hpp"

#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "sequence/seq.hpp"

using namespace std;

class BadConversion : public std::runtime_error {
  public:
	BadConversion(std::string const& s) : std::runtime_error(s) { }
};

inline double convertToInt(std::string const& s) {
	istringstream i(s);
	int x;
	if (!(i >> x))
		throw BadConversion("convertToInt(\"" + s + "\")");
	return x;
}

int main(int argc, char * argv[]) {
	if (argc < 5) {
		cout << "Usage: ./reconstruct ufFilename readFilename nthreads qvoffset [outFilename=\"reads.processed\"]\n";
		return 0;
	}
	cout << argc << "\n";
	string ufFilename = argv[1];
	string readsFilename = argv[2];
	int nthreads = atoi(argv[3]);
	int qvoffset = atoi(argv[4]);
	string outFilename = "reads.processed";
	if (argc > 5) outFilename = argv[5];
	cout << "Starting. uf = " << ufFilename << ". Threads = " << nthreads << ".\n";
	flush(cout);

	ifstream inf;
	open_file(inf, ufFilename);
	// map of kmers
	map<KMer, KMer, KMer::less2> changes;
	int blockNum = 0;
	int curBlockNum = 0;
	KMer curCenter; bool hasCenter = false;
	vector<KMer> block;
	vector<string> row;
	size_t count = 0;
	while (get_row_whitespace(inf, row)) {
		++count; if (count % 1000000 == 0) cout << count << "\n";
		if (row.size() < 1) continue;
		if (row[6] == "goodSingleton") continue;
		curBlockNum = convertToInt(row[0]);
		
		if (curBlockNum != blockNum) { // a block is over
			if (!hasCenter && block.size() > 0) {
				cout << "Block " << blockNum << " without center!!!\n";
			} else {
				for (size_t i=0; i<block.size(); ++i) {
					changes.insert(make_pair(block[i], curCenter));
				}
				blockNum = curBlockNum; block.clear(); hasCenter = false;
			}
		}
		if (row[6] == "center") {
			curCenter = KMer(row[1]); hasCenter = true;
		} else if (row[6] == "change") {
			block.push_back(KMer(row[1]));
		}
	}
	inf.close();

	// now change the reads
	ireadstream ifs(readsFilename.data(), qvoffset);
	ofstream ofs;
	ofs.open(outFilename);
	Read r;
	map<KMer, KMer, KMer::less2>::iterator it;
	while (!ifs.eof()) {
		ifs >> r;
		// trim the reads for bad quality and process only the ones with at least K "reasonable" elements
		if (r.trimBadQuality() >= K) {
			string seq = r.getSequenceString();
			// create auxiliary structures for consensus
			vector<int> vA(r.size(), 0), vC(r.size(), 0), vG(r.size(), 0), vT(r.size(), 0);
			vector< vector<int> > v;  // A=0, C=1, G=2, T=3
			v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
			
			size_t i = 0;
			bool changedRead = false;
			while (true) {
				i = r.firstValidKmer(i, K);
				if (i+K > r.size()) break;
				KMer kmer = KMer(r.getSubSequence(i, K));
				it = changes.find(kmer);
				
				if (it != changes.end()) { // we've found that we have to change this kmer
					// pretty print the k-mer
					if (!changedRead) {
						changedRead = true;
						cout << "\n " << r.getName() << "\n" << r.getSequenceString().data() << "\n";
					}
					for (size_t j=0; j<i; ++j) cout << " ";
					cout << it->second.str().data() << "\n";

					for (size_t j=0; j<K; ++j) {
						if (kmer[j] != it->second[j]) {
							v[it->second[j]][i+j]++;
						}
					}
				}
				
				// go to next kmer
				if (i+K < r.size() && is_nucl(seq[i+K])) {
					kmer = kmer << r[i+K];
					++i;
				} else {
					i = i+K;
					break;
				}
			}
				
			// at this point the array v contains votes for consensus
			
			
			// find max consensus element
			for (size_t j=0; j<r.size(); ++j) {
				char cmax = seq[j]; int nummax = 0;
				for (size_t k=0; k<4; ++k) {
					if (v[k][j] > nummax) {
						cmax = nucl(k); nummax = v[k][j];
					}
				}
				seq[j] = cmax;
			}
			if (changedRead) {
 				for (size_t i=0; i<4; ++i) {
					for (size_t j=0; j<r.size(); ++j) {
						if (v[i][j] > 9) cout << "*"; else cout << v[i][j];
					}
					cout << "\n";
				}
				cout << seq.data() << "\n";
			}


			ofs << "@" << r.getName() << endl << seq.data() << endl << "+" << r.getName() << endl << r.getPhredQualityString(qvoffset) << endl;
		}
	}
	ofs.close();
	inf.close();
	
	return 0;
}


/*
 * correct.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#include <unordered_map>

#include "defs.hpp"
#include "read.hpp"
#include "logging.hpp"
#include "mathfunctions.hpp"
#include "ireadstream.hpp"
#include "uf.hpp"
#include "common/sequence.hpp"
#include "hammer_config.hpp"

using namespace std;
using namespace __gnu_cxx;


LOGGER("h");

/**
 * correct one read
 * can there be cycles and other dragons? who knows...
 */
string correctSequence(const Sequence &s, const KMerHashMap & m) {
	string res = s.str();
	cout << res << "  " << res.size() << "\n";
	if (res.size() < K) return res;
	bool changed = true;
	int start = 0;
	while (changed) {
		changed = false;
		KMer kmer(res.substr(start, start+K));
		for (size_t i=start; i<=res.size()-K; ++i) {
			cout << "looking for " << kmer << endl;
			KMerHashMap::const_iterator it = m.find(kmer);
			if (it != m.end()) {
				cout << "  found " << it->second.str().data() << "\n";
				changed = true;
				for (size_t j=0; j<K; ++j)
					res[i+j] = nucl(it->second[j]);
				start = max(0, (int)i-K+1);
				break;
			}
			kmer = kmer << res[i+K];
		}
	}
	return res;
}


int main(int argc, char * argv[]) {
	string ufFilename = argv[1];
	string readsFilename = argv[2];
	string outFilename = argv[3];
	
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

	ireadstream ifs(readsFilename.data());
	ofstream ofs;
	ofs.open(outFilename.data());
	Read r;
	while (!ifs.eof()) {
		ifs >> r;
		string corrs = correctSequence(r.getSequence(), ufc_map);
		ofs << "@" << r.getName() << endl << corrs.data() << endl << "+" << endl << r.getPhredQualityString() << endl;
	}

	ifs.close(); ofs.close();
	return 0;
}


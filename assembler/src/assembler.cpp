#include <cassert>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include "parser.hpp"
//#include "condensed_graph.hpp"

using namespace std;

pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");

#define MPSIZE 100

int main(int argc, char *argv[]) {
  
	std::cerr << "Hello, I am assembler!" << std::endl;
	time_t now = time(NULL);
	//SeqVarLen x = SeqVarLen("ACGTACGTTCGATCGAAA");
	//cerr << x.str() << endl;
	//return 0;
	FASTQParser<MPSIZE>* fqp = new FASTQParser<MPSIZE>();
	fqp->open(filenames.first, filenames.second);
	vector<MatePair<MPSIZE> > mps;
	int cnt = 0;
	while (!fqp->eof()) {
		MatePair<MPSIZE> mp = fqp->read(); // is it copy? :)
		if (mp.id != -1) { // have 'N' in reads
			mps.push_back(mp); // it's copy! :)
			//cout <<  mp.id << endl << mp.seq1.str() << endl <<  mp.seq2.str() << endl;
		}
		cnt++;
	}
	cout << "Total reads: " << cnt << endl;
	cout << "Clear (without N) reads: " << mps.size() << endl;
	cout << "seconds: " << (time(NULL) - now) << endl;
	fqp->close();
	return 0;
}

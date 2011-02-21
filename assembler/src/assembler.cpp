#include <cassert>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include "parser.hpp"

using namespace std;

pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");

int main(int argc, char *argv[]) {
  
	std::cerr << "Hello, I am assembler!" << std::endl;
	time_t now = time(NULL);

	FASTQParser* fqp = new FASTQParser();
	fqp->open(filenames.first, filenames.second);
	vector<MatePair> mps;
	while (!fqp->eof()) {
		MatePair mp = fqp->read(); // is it copy? :)
		if (mp.id != -1) { // have 'N' in reads
			mps.push_back(mp); // it's copy! :)
			cout <<  mp.id << endl << mp.seq1.str() << endl <<  mp.seq2.str() << endl;
		}
	}
	cout << mps.back().id << endl;
	cout << mps.size() << endl;
	cout << "seconds: " << (time(NULL) - now) << endl;
	fqp->close();
	return 0;
}

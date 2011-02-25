#include <cassert>
#include <iostream>
#include <list>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include "../parser.hpp"
//#include "condensed_graph.hpp"

using namespace std;

pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");

#define MPSIZE 100

int main(int argc, char *argv[]) {

	std::cerr << "Hello, I am assembler!" << std::endl;
	time_t now = time(NULL);
	Seq<10> x = Seq<10>("ACGTACGTTC");
	Seq<10> y = x.shift_right(0);
	Seq<10> z = y.shift_left(3);
	cerr << "Seq Shift Test 1: " << x.str() << endl;
	cerr << "Seq Shift Test 2: " << y.str() << endl;
	cerr << "Seq Shift Test 3: " << z.str() << endl;
	cerr << "sizeof(Seq<10>): " << sizeof(z) << endl;
	return 0;
	FASTQParser<MPSIZE>* fqp = new FASTQParser<MPSIZE>();
	fqp->open(filenames.first, filenames.second);
	vector<MatePair<MPSIZE> > mps;
	int cnt = 0;
	while (!fqp->eof()) {
		MatePair<MPSIZE> mp = fqp->read(); // is it copy? :)
		if (mp.id != -1) { // don't have 'N' in reads
			mps.push_back(mp); // it's copy! :)
			//cout <<  mp.id << endl << mp.seq1.str() << endl <<  mp.seq2.str() << endl;
		}
		cnt++;
	}
	cout << "Total reads: " << cnt << endl;
	cout << "Clear (without N) reads: " << mps.size() << endl;
	cout << "seconds: " << (time(NULL) - now) << endl;
	fqp->close();
	freopen("./data/reads.txt", "w", stdout);
	for(int i = 0; i < mps.size(); i++) {
		Seq<MPSIZE> tmp = mps[i].seq1;
		string s = tmp.str();
		std::cout<<s;
		tmp = mps[i].seq2;
		s = tmp.str();
		std::cout<<" "<<s<<endl;
		//std::cout << mps[i].seq1::str() << mps[i].seq2::str() << endl;
	}
	return 0;
}

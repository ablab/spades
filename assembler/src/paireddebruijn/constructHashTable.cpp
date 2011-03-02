#include "common.hpp"
#include "../sequence.hpp"
#include "constructHashTable.hpp"

using namespace std;

typedef map<ll, vector<ll> > myMap;

typedef vector<Sequence*> downSeqs;


int totalKmers = 0;
int uniqPairs = 0;
ll upperMask = (((ll) 1) << (2 * k)) - 1;
ll lowerMask = (((ll) 1) << (2 * l)) - 1;

ll upperMax = ((ll) 1) << 46;
const int MAXLMERSIZE = 10000;

string decompress(ll a, int l) {

	string res = "";
	res.reserve(l);
	forn(i,l)
		res += " ";
	forn(i, l) {
		res[l - i - 1] = nucl((a & 3));
		a >>= 2;
	}
	return res;
}

inline int codeNucleotide(char a) {
	if (a == 'A')
		return 0;
	else if (a == 'C')
		return 1;
	else if (a == 'G')
		return 2;
	else if (a == 'T')
		return 3;
	else {
		std::cerr << "oops!";
		return -1;
	}
}
/*void testSequence(){
	srand(239);
	forn(i, 1000) {

		ll ts = ((ll) rand()) * ((ll) rand());
		//cerr << ts;
		string s = decompress(ts, l);

		Sequence* tst = new Sequence(s);
		string ss = tst->Str();
		assert (ss == s);
		//cerr << s <<endl<< ss<<endl<<endl;

	}
}*/
void codeRead(char *read, char *code) {
	for (int i = 0; i < readLength; i++) {
		code[i] = codeNucleotide(read[i]);
	}
}

//toDo
downSeqs clusterize(ll* a, int size, int max_shift) {
	downSeqs res;
	res.clear();
	assert (max_shift <= 20);
	int right[MAXLMERSIZE];
	int left[MAXLMERSIZE];
	int used[MAXLMERSIZE];
	int shift_left[MAXLMERSIZE];
	int shift_right[MAXLMERSIZE];
	//-1 = no neighbor;
	//-2 = more than 1 neighbor
	forn(i, size) {
		right[i] = -1;
		left[i] = -1;
		used[i] = 0;
		shift_left[i] = 0;
		shift_right[i] = 0;
	}
	ll diff;
	forn(i, size) {
		ll right_tmp = a[i];
		ll left_tmp = a[i];
		ll p2 = 0;
		ll upper_bound;
		forn(shift, max_shift) {
		    right_tmp = ((right_tmp << 2) & lowerMask);
		    p2 += 2;
		    upper_bound = ((ll) 1) << p2;
		    if (!shift_right[i]){
				forn(j, size) {
					diff = a[j] - right_tmp;
					if ((diff >= 0) && (diff < upper_bound) && (i != j)){
						shift_right[i] = p2/2;
						if (right[i] == -1) {
							right[i] = j;
						}
						else
							right[i]  = -2;
					}
				}
		    }
			left_tmp >>= 2;
			if (!shift_left[i]) {
				forn(j, size) {
					diff = a[j] - left_tmp;
					if ((i != j) && ((diff & (lowerMask >> p2)) == 0)){
						shift_left[i] = p2/2;
						if (left[i] == -1)
							left[i] = j;
						else
							left[i]  = -2;
					}
				}
			}
		}
	}
	int color = 1;
	forn(i, size) {
		int seqlength = l;
		if (used[i] == 0) {
			int ii = i;
			used[i] = color;
			while ((left[ii] >= 0) && (right[left[ii]] == ii) && (left[ii] != i)){
				seqlength += shift_left[ii];
				ii = left[ii];
				used[ii] = color;
			}
			int leftend = ii;

			ii = i;
			while ((right[ii] >= 0) && (left[right[ii]] == ii) && (right[ii] != i)){
				seqlength += shift_right[ii];
				ii = right[ii];
				used[ii] = color;
				seqlength++;
			}
			int rightend = ii;
			ii = leftend;
			string s = decompress(a[leftend], l);
			while (ii != rightend) {
		//		cerr << "clusterizing....";
				int p = shift_right[ii];
				ll maxsd = ((ll) 3) << (2 * (p-1));
				ii = right[ii];
				forn(j, p)
					s += nucl((a[ii] & maxsd) >> (2*(p-j-1)));
			}
			Sequence* tmpSeq = new Sequence(s);
			res.pb(tmpSeq);
			color++;
		}
	}
	/*{
		forn(i, size) {
			cerr << left[i] << " ";
		}
		cerr << endl;
		forn(i, size) {
			cerr << right[i] << " ";
		}
	}*/
//	assert(0);
	return res;
}

inline bool nextReadPair(char * &read1, char * &read2) {
	return scanf("%s %s", read1, read2) == 2;
}

ll extractMer(char *read, int shift, int length) {
	ll res = 0;
	for (int i = 0; i < length; i++) {
		res = res << 2;
		res += read[shift + i];
	}
	return res;
}

inline bool checkBoundsForUpper(ll upper) {
	return true;
	if ((upper >= 1<<20) && (upper < upperMax))
		return true;
	else return false;
}

void addPairToTable(myMap& table, ll upper, ll lower) {
	if (table.find(upper) != table.end()) {
		if (find(table[upper].begin(), table[upper].end(), lower)
				== table[upper].end()) {
			table[upper].pb(lower);
			++uniqPairs;
		}
	} else {
		vector<ll> tmp;
		tmp.clear();
		tmp.pb(lower);
		table.insert(make_pair(upper, tmp));

//		cerr<<"inserting"<<table.size();
		++uniqPairs;
	}
}

void processReadPair(myMap& table, char *upperRead, char *lowerRead) {
	int shift = (l - k) / 2;
	ll upper = extractMer(upperRead, shift, k);
	ll lower = extractMer(lowerRead, 0, l);
//	fprintf(stderr,"%lld %lld\n", upper, lower);
	for (int j = 0; j + l <= readLength; j++) {
		if (checkBoundsForUpper(upper)) {
			addPairToTable(table, upper, lower);
			totalKmers++;

		}

		upper <<= 2;
		upper += upperRead[j + l - shift];
		upper &= upperMask;

		lower <<= 2;
		lower += lowerRead[j + l];
		lower &= lowerMask;
		//fprintf(stderr,"%d %d\n", upper, lower);

	}
//	cerr << table.size()<<endl;
}


void constructTable(myMap &table) {
	int count = 0;
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	while (nextReadPair(upperNuclRead, lowerNuclRead)) {
//		fprintf(stderr, "%s", upperNuclRead);
		codeRead(upperNuclRead, upperRead);
		codeRead(lowerNuclRead, lowerRead);
		processReadPair(table, upperRead, lowerRead);
		if (!(count & (1024*128 - 1)))
			cerr<<"read number "<<count<<" processed"<<endl;
		count++;
	}
}

void outputTable(myMap &pairedTable) {
	int j = 0;
	for (myMap::iterator iter = pairedTable.begin() ; iter != pairedTable.end(); iter++) {
		pair<ll, vector<ll> > p = (*iter);
		cout << p.fi << " " << p.se.size() << endl;
		forn(i, p.se.size()) {
			cout << p.se[i] << " ";
		}
		cout << endl << endl;
		if (!(j & 1023))
			cerr << j << endl;
		j++;
	}
	pairedTable.clear();
}

void readsToPairs(string inputFile, string outputFile) {

	myMap table;
	cerr << "generation of k-l pairs started"<<endl;
	freopen(inputFile.c_str(), "r", stdin);
	constructTable(table);
	cerr << "generation of k-l pairs finished, dumping to disk."<<endl;
	freopen(outputFile.c_str(), "w", stdout);
	cerr<< "outputFile opened";
	outputTable(table);
	fclose(stdout);
	table.clear();
}
//#define OUTPUT_DECOMPRESSED
int pairsToSequences(string inputFile, string outputFile) {
	FILE* inFile = freopen(inputFile.c_str(), "r", stdin);
#ifdef OUTPUT_DECOMPRESSED
	FILE* decompressed = fopen("data/decompressed.out", "w" );
#endif
	int ok = 1;
	ll lmers[MAXLMERSIZE];
	ll kmer;
	int lsize;
	FILE* outFile = fopen(outputFile.c_str(), "w");
	int count = 0;
	while (1) {
		count++;
		ok = fscanf(inFile, "%lld %d", &kmer, &lsize);
		if (ok != 2) {
			if (ok != 0)
				cerr<< "error in reads.";
			else
				cerr << "Finished!!";
			break;
		}
		if (lsize > MAXLMERSIZE) {
			cerr << "TOO BIIIIG";
			return -2;
		}

		forn(i, lsize) {

			if (fscanf(inFile, "%lld", &lmers[i]) != 1) {
				cerr << "Error in main1 reading l-mers";
				return -1;
			}
		}
		sort(lmers, lmers + lsize);
		downSeqs clusters =  clusterize(lmers, lsize, 3);
		int clsize = clusters.size();
		string outstring;

//		string s = decompress(kmer, k);
//		fprintf(decompressed, "%s %d\n", s.c_str(), lsize);
		fprintf(outFile, "%lld %d\n", kmer, clsize);
#ifdef OUTPUT_DECOMPRESSED
		forn(i, lsize) {
			fprintf(decompressed, "%s ", decompress(lmers[i], l).c_str());
		}
#endif
		forn(i, clsize) {
			outstring = clusters[i]->Str();
			fprintf(outFile, "%s ",outstring.c_str());
		}
		fprintf(outFile, "\n");
#ifdef OUTPUT_DECOMPRESSED
		fprintf(decompressed, "\n");
#endif
		//	return 0;
		if (!(count & ((1 << 15) - 1) ))
			cerr<< "klmer numero "<< count <<"generated" <<endl;
	}
	cerr<<"finished";
	return 0;
}

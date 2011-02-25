#include "common.h"
#include "../seq.hpp"
#include "constructHashTable.hpp"

using namespace std;

typedef map<ll, vector<ll> > myMap;

typedef vector<Sequence> downSeqs;
int k = 25;
int l = 31;
int readLength = 100;
long totalKmers = 0;
long uniqPairs = 0;
ll upper_max = ((ll) 1) << 46;

inline int codeNucleotide(char a) {
	if (a == 'A')
		return 0;
	else if (a == 'C')
		return 1;
	else if (a == 'T')
		return 2;
	else if (a == 'G')
		return 3;
	else {
		std::cerr << "oops!";
		return -1;
	}
}

void codeRead(char *read, char *code) {
	for (int i = 0; i < readLength; i++) {
		code[i] = codeNucleotide(read[i]);
	}
}

string decompress(ll a) {
	string res = "";
	res.reserve(l);
	forn(i,l)
		res += " ";
	forn(i, l) {
		res[l - i - 1] = '0' + (a & 3);
		a >>= 2;
	}
	return res;
}

//toDo
downSeqs clusterize(int* a, int size) {
	downSeqs res;
	return res;
}

bool nextReadPair(char * &read1, char * &read2) {
	return scanf("%s %s", read1, read2) == 2;
}

ll extractMer(char *read, int shift, int length) {
	ll res = 0;
	for (int i = 0; i < length; i++) {
		res = res << 2;
		res += read[shift + length];
	}
	return res;
}

bool checkBoundsForUpper(ll upper) {
	return true;
}

void addPairToTable(myMap table, ll upper, ll lower) {
	if (table.find(upper) != table.end()) {
		if (find(table[upper].begin(), table[upper].end(), lower)
				== table[upper].end()) {
			table[upper].pb(lower);
			++uniqPairs;
		}
	} else {
		vector<ll> tmp;
		tmp.pb(lower);
		table.insert(make_pair(upper, tmp));
		++uniqPairs;
	}
}

void processReadPair(myMap table, char *upperRead, char *lowerRead) {
	ll upperMask = (((ll) 1) << (2 * k)) - 1;
	ll lowerMask = (((ll) 1) << (2 * l)) - 1;
	int shift = (l - k) / 2;
	ll upper = extractMer(upperRead, shift, k);
	ll lower = extractMer(lowerRead, 0, l);
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
	}
}

void constructTable(myMap &table) {
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	while (nextReadPair(upperNuclRead, lowerNuclRead)) {
		codeRead(upperNuclRead, upperRead);
		codeRead(lowerNuclRead, lowerRead);
		processReadPair(table, upperRead, lowerRead);
	}
}

void outputTable(myMap pairedTable) {
	int j = 0;
	for (myMap::iterator iter = pairedTable.begin(); iter != pairedTable.end(); iter++) {
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

void readsToPairs(char *inputFile, char *outputFile) {
	myMap table;
	freopen(inputFile, "r", stdin);
	constructTable(table);
	freopen(outputFile, "w", stdout);
	outputTable(table);
	fclose(stdout);
	delete &table;
}

int main1() {
	FILE* f = freopen("data/reads.out", "r", stdin);
	cerr << f << endl;
	int ok = 1;
	const int MAXLMERSIZE = 2000;
	ll lmers[MAXLMERSIZE];

	ll kmer;
	int lsize;
	while (1) {
		ok = scanf("%lld %d", &kmer, &lsize);
		if (ok != 2) {
			break;
		}
		forn(i, lsize) {
			scanf("%lld", &lmers[i]);
			//	cerr << lmers[i]<<endl;
		}
		sort(lmers, lmers + lsize);
		string s = decompress(lmers[0]);
		forn(i, lsize)
			cerr << lmers[i] << ":" << decompress(lmers[i]) << " ";
		cerr << endl << endl;
	}
}

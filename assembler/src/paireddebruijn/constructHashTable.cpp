#include "common.h"

using namespace std;

//typedef __gnu_cxx::hash_map<ll, vector<ll> > myMap;
typedef map<ll, vector<ll> > myMap;

myMap pairedTable;
int k;
int l;
int read_length;

inline int value(char a) {
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
int main() {
	//freopen("config.ini", "r", stdin);
	//scanf ("Upper k-mer size = %d",&k);
	//scanf ("Lower k-mer size = %d",&l);
	freopen("reads.txt", "r", stdin);
	l = 31;
	k = 25;
	ll upper_cut = (((ll) 1) << (2 * k)) - 1;

	ll lower_cut = (((ll) 1) << (2 * l)) - 1;
	read_length = 100;
	int shift = (l - k) / 2;
	int maxn = 1 << 20;
	int read_num = 0;
    
	long totalKmers=0;
	long uniqPairs=0;
	ll upper_max = ((ll) 1) << 55;

	while (1) {
		if (!(read_num & 1023))
			cerr << "read:" << read_num <<"  Lmers: "<<totalKmers<<  "Unique: "<<uniqPairs<<endl;
		read_num++;
		char r1[102];
		char r2[102];
		char rr1[102];
		char rr2[102];
		int n = scanf("%s %s", &rr1, &rr2);
		if (n != 2)
			break;
		//cerr << n;
		//cerr<< rr1;
		//return 0;
		forn(i, read_length) {
			r1[i] = value(rr1[i]);
			r2[i] = value(rr2[i]);
		}
		ll upper = 0;
		ll lower = 0;
		forn(j, k) {
			upper = upper << 2;
			upper += r1[j + shift];
		}
		forn(j, l) {
			lower = lower << 2;
			lower += r2[j];
		}
		forn(j, read_length - l+1) {
//			if ((upper > 0) && (upper < upper_max)) {
			if (1){
				if (pairedTable.find(upper) != pairedTable.end()) {
					if(find(pairedTable[upper].begin(), pairedTable[upper].end(), lower) == pairedTable[upper].end())
					    {pairedTable[upper].pb(lower);++uniqPairs;}
				}
				else {
					vector<ll> tmp;
					tmp.pb(lower);
					pairedTable.insert(make_pair(upper, tmp));
					++uniqPairs;
				}
			totalKmers++;
			}
			upper <<= 2;
			lower <<= 2;

			upper += r1[j + l - shift];
			lower += r2[j + l];

			upper &= upper_cut;
			lower &= lower_cut;
		}
	}
	freopen("./data/reads.out", "w", stdout);
	int j = 0;
	for (myMap::iterator iter = pairedTable.begin(); iter != pairedTable.end(); iter++) {
		pair<ll, vector<ll> > p = (*iter);
		//		cerr<<j<<endl;
		cout << p.fi << " "<< p.se.size() <<endl;
		//		cerr << p.fi << endl;
		forn(i, p.se.size()) {
			cout << p.se[i] << " ";
		}
		cout << endl << endl;
		if (!(j & 1023))
			cerr << j << endl;
		j++;
	}
}

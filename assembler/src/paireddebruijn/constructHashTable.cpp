#include "common.h"


using namespace std;

typedef __gnu_cxx::hash_map<int, vector<int> > myMap;
myMap pairedTable;
int k;
int l;
int read_length;

inline int value(char a){
	if (a == 'A') return 0;
	else if (a == 'C') return 1;
	else if (a == 'T') return 2;
	else if (a == 'G') return 3;
	else {
		std::cerr << "oops!";
		exit;
	}

}
int main(){
	freopen("config.ini", "r", stdin);
	scanf ("Upper k-mer size = %d",&k);
	scanf ("Lower k-mer size = %d",&l);
	freopen("reads.txt", "r", stdin);
	l = 31;
	k = 25;
	ll upper_cut = (((ll) 1) << (2 * l)) - 1;

	ll lower_cut = (((ll) 1) << (2 * k)) - 1;
	read_length = 100;
	int shift = (l - k) / 2;
	int maxn = 1000000;
	forn(i, maxn) {

		char r1[102];
		char r2[102];
		char rr1[102];
		char rr2[102];
		scanf("%s %s", rr1, rr2);
		forn(i, read_length) {
			r1[i] = value(rr1[i]);
			r2[i] = value(rr2[i]);
		}
		ll upper = 0;
		ll lower = 0;
		forn(j, l) {
			upper = upper << 2;
			upper += r1[j];
		}
		forn(j, k) {
			lower = lower << 2;
			lower += r1[j + shift];
		}
		forn(j, read_length - l) {
			if (pairedTable.find(upper) != pairedTable.end())
				pairedTable[upper].pb(lower);
			else {
				vector<ll> tmp;
				tmp.pb(lower);
				pairedTable.insert(make_pair(lower,tmp));
			}
			upper <<=2;
			lower <<=2;

			upper += r1[j + l - shift];
			lower += r2[j + l];

			upper &= upper_cut;
			lower &= lower_cut;
		}
	}
	for (myMap::iterator iter = pairedTable.begin(); iter !=pairedTable.end(); iter++ ) {
		pair<ll, vector<ll> > p = (*iter);
		cerr << p.fi << endl;
		forn(i, p.se.size()){
			cerr<< p.se[i]<<" ";
		}

	}
}

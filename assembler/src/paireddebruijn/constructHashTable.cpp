#include "common.hpp"
#include "sequence.hpp"
#include "constructHashTable.hpp"
#include "graphio.hpp"
LOGGER("p.constructHashTable");

using namespace std;

typedef vector<pair<Sequence*,int>> downSeqs;


int totalKmers = 0;
int uniqPairs = 0;
const int MAXLMERSIZE = 10000;
ll upperMask;
ll lowerMask;

ll upperMax;


//toDo
void initGlobal(){
	upperMask = (((ll) 1) << (2 * k)) - 1;
	lowerMask = (((ll) 1) << (2 * l)) - 1;
	upperMax = ((ll) 1) << 46;
}

/*
downSeqs oldclusterizeLset(pair<ll,int>* a, int size, int max_shift, set<ll> &lset) {
	downSeqs res;
	res.clear();
	assert (max_shift <= 20);
//	cerr << lset.size()<<endl;
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
		ll right_tmp = a[i].first;
		ll left_tmp = a[i].first;
		ll p2 = 0;
		ll upper_bound;
		forn(shift, max_shift) {
		    right_tmp = ((right_tmp << 2) & lowerMask);
		    p2 += 2;
		    int cright = 0;
		  /*  if (!shift_right[i]) {
		    	forn(ii, (1<<p2)) {
		    		if (lset.find(ii + right_tmp) != lset.end()) {
		    			cright ++;
		    		}
		    		if (cright > 1) {
		    			shift_right[i] = -2;
		    			break;
		    		}
		    	}
			}

		    //cerr <<"cright" <<cright << endl;

		    if (!(cright == 0 || shift_right[i] || cright > 1)) {


		    	upper_bound = ((ll) 1) << p2;
				forn(j, size) {
					diff = a[j].first - right_tmp;
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
			cright = 0;
			if (!shift_left[i]) {
				forn(ii, (1<<p2)) {
					ll left_n = ((ll) ii) << (2*l - p2);
	//				cerr<<a[i] << " "<< left_n + left_tmp <<endl;
					if (lset.find(left_n + left_tmp) != lset.end()) {
						cright ++;
					}
					if (cright > 1) {
						shift_left[i] = -2;
						break;
					}
				}
			}
		//	cerr <<"cleft" <<cright << endl;
			if (!(cright == 0 || shift_left[i] || cright > 1)) {
				forn(j, size) {
					diff = a[j].first - left_tmp;
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
			//cerr <<"color = :"<< color << endl;
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
			string s = decompress(a[leftend].first, l);
			int coverage = a[leftend].second;
			while (ii != rightend) {
		//		cerr << "clusterizing....";
				int p = shift_right[ii];
				ll maxsd = ((ll) 3) << (2 * (p-1));
				ii = right[ii];
				forn(j, p) {
				//	cerr << ((a[ii] & maxsd) >> (2*(p-j-1)));
					s += nucl((a[ii].first & maxsd) >> (2*(p-j-1)));
					maxsd >>= 2;
			//		cerr << "OK" <<endl;
					if (coverage < a[ii].second) coverage = a[ii].second;
				}
			}
			Sequence* tmpSeq = new Sequence(s);
			res.pb(make_pair(tmpSeq,coverage));
			color++;
		}
	}
	/*
	{
		forn(i, size) {
			cerr << left[i] << " ";
			cerr << shift_left[i] << " ";

		}
		cerr << endl;
		forn(i, size) {
			cerr << right[i] << " ";
			cerr << shift_right[i] << " ";
		}
		forn(i, res.size()) {
			cerr<<(res[i].first)->str() << endl;
		}
	}
	assert(0);*//*
	return res;
}
*/



downSeqs clusterize(pair<ll,int>* a, int size, int max_shift) {
	downSeqs res;
	res.clear();
	vector<string> tmp_res;

	vector<int> tmp_cov;
	assert (max_shift <= 20);
//	cerr << lset.size()<<endl;
	int right[MAXLMERSIZE];
	int left[MAXLMERSIZE];
	int used[MAXLMERSIZE];
	int shift_left[MAXLMERSIZE];
	int shift_right[MAXLMERSIZE];
	//-1 = no neighbor;
	//-2 = more than 1 neighbor
	DEBUG("clusterizing");
	forn(i, size) {
		DEBUG(decompress(a[i].first, l)<< " "<< i);
		right[i] = -1;
		left[i] = -1;
		used[i] = 0;
		shift_left[i] = 0;
		shift_right[i] = 0;
	}
	ll diff;
	forn(i, size) {
		ll right_tmp = a[i].first;
		ll left_tmp = a[i].first;
		ll p2 = 0;
		if (a[i].second < 1) {
//			used[i] = true;
//			continue;
		}
		ll upper_bound;
		forn(shift, max_shift) {
		    right_tmp = ((right_tmp << 2) & lowerMask);
		    p2 += 2;
		    int cright = 0;
		    if (!( shift_right[i] )) {
		    	upper_bound = ((ll) 1) << p2;
				forn(j, size) {
					diff = a[j].first - right_tmp;
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
			cright = 0;
			if ( !(shift_left[i] )) {
				forn(j, size) {
					diff = a[j].first - left_tmp;
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
	vector<int> leftway;
	forn(i, size) {
		int seqlength = l;
		if (used[i] == 0) {
			int ii = i;
			leftway.clear();
			DEBUG("COLOR: " << color << " from i: "<< i);
			used[i] = color;
			//cerr <<"color = :"<< color << endl;
			while ((left[ii] >= 0) && (left[ii] != i)){
				seqlength += shift_left[ii];
				leftway.pb(ii);
				ii = left[ii];
				used[ii] = color;
				DEBUG(ii);
			}
			int leftend = ii;
				DEBUG("righrt");
			ii = i;
			while ((right[ii] >= 0) && (right[ii] != i)){
				seqlength += shift_right[ii];
				ii = right[ii];
				used[ii] = color;
				seqlength++;
				DEBUG(ii);
			}
			int rightend = ii;
			ii = leftend;
			string s = decompress(a[leftend].first, l);
			int coverage = a[leftend].second;
			DEBUG("leftstirng "<< s);
			int currInd = leftway.size()-1;
			while (ii != rightend) {
		//		cerr << "clusterizing....";
				int p = shift_right[ii];
				ll maxsd = ((ll) 3) << (2 * (p-1));
				if (currInd >= 0) {
					ii = leftway[currInd];
					currInd --;
				}
				else
					ii = right[ii];
				DEBUG(ii << " " << p <<" " << maxsd);
				forn(j, p) {
				//	cerr << ((a[ii] & maxsd) >> (2*(p-j-1)));
					s += nucl((a[ii].first & maxsd) >> (2*(p-j-1)));
					maxsd >>= 2;
			//		cerr << "OK" <<endl;
					if (coverage < a[ii].second) coverage = a[ii].second;
				}
			}
			DEBUG("seq: s" << s);
			tmp_res.push_back(s);
			tmp_cov.pb(coverage);
//			Sequence* tmpSeq = new Sequence(s);
//			res.pb(make_pair(tmpSeq,coverage));
			color++;
		}
	}
	forn(i,  tmp_res.size()) {
		int good = 1;
		for(int j = 0; j<tmp_res.size(); j++){
			if (j != i && tmp_res[j].find(tmp_res[i]) != string::npos) {
				good = 0;
				break;
				INFO("SUBSEQ" << tmp_res[i] << " " << tmp_res[j]);
			}
		}
		if (good) {
			Sequence* tmpSeq = new Sequence(tmp_res[i]);
			res.pb(make_pair(tmpSeq,tmp_cov[i]));

		}
	}
	/*
	{
		forn(i, size) {
			cerr << left[i] << " ";
			cerr << shift_left[i] << " ";

		}
		cerr << endl;
		forn(i, size) {
			cerr << right[i] << " ";
			cerr << shift_right[i] << " ";
		}
		forn(i, res.size()) {
			cerr<<(res[i].first)->str() << endl;
		}
	}
	assert(0);*/
	return res;
}


downSeqs clusterize0704(pair<ll,int>* a, int size, int max_shift) {
	downSeqs res;
	res.clear();
	vector<string> tmp_res;

	vector<int> tmp_cov;
	assert (max_shift <= 20);
//	cerr << lset.size()<<endl;
	int right[MAXLMERSIZE];
	int left[MAXLMERSIZE];
	int used[MAXLMERSIZE];
	int shift_left[MAXLMERSIZE];
	int shift_right[MAXLMERSIZE];
	//-1 = no neighbor;
	//-2 = more than 1 neighbor
	DEBUG("clusterizing");
	forn(i, size) {
		DEBUG(decompress(a[i].first, l)<< " "<< i);
		right[i] = -1;
		left[i] = -1;
		used[i] = 0;
		shift_left[i] = 0;
		shift_right[i] = 0;
	}
	ll diff;
	forn(i, size) {
		ll right_tmp = a[i].first;
		ll left_tmp = a[i].first;
		ll p2 = 0;
		if (a[i].second < 1) {
//			used[i] = true;
//			continue;
		}
		ll upper_bound;
		forn(shift, max_shift) {
		    right_tmp = ((right_tmp << 2) & lowerMask);
		    p2 += 2;
		    int cright = 0;
		    if (!( shift_right[i] )) {
		    	upper_bound = ((ll) 1) << p2;
				forn(j, size) {
					diff = a[j].first - right_tmp;
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
			cright = 0;
			if ( !(shift_left[i] )) {
				forn(j, size) {
					diff = a[j].first - left_tmp;
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
	vector<int> leftway;
//	memset()
	forn(i, size) {

		int seqlength = l;
		if (used[i] == 0) {
			int ii = i;
			int intersect = -1;
			leftway.clear();
			DEBUG("COLOR: " << color << " from i: "<< i);
			used[i] = color;
			//cerr <<"color = :"<< color << endl;
			while ((left[ii] >= 0) && (left[ii] != i)){
				seqlength += shift_left[ii];
				leftway.pb(ii);
				ii = left[ii];
				used[ii] = color;
				DEBUG(ii);
			}
			int leftend = ii;
				DEBUG("righrt");
			ii = i;
			while ((right[ii] >= 0) && (right[ii] != i)){
				seqlength += shift_right[ii];
				ii = right[ii];
				used[ii] = color;
				seqlength++;
				DEBUG(ii);
			}
			int rightend = ii;
			ii = leftend;
			string s = decompress(a[leftend].first, l);
			int coverage = a[leftend].second;
			DEBUG("leftstirng "<< s);
			int currInd = leftway.size()-1;
			while (ii != rightend) {
		//		cerr << "clusterizing....";
				int p = shift_right[ii];
				ll maxsd = ((ll) 3) << (2 * (p-1));
				if (currInd >= 0) {
					ii = leftway[currInd];
					currInd --;
				}
				else
					ii = right[ii];
				DEBUG(ii << " " << p <<" " << maxsd);
				forn(j, p) {
				//	cerr << ((a[ii] & maxsd) >> (2*(p-j-1)));
					s += nucl((a[ii].first & maxsd) >> (2*(p-j-1)));
					maxsd >>= 2;
			//		cerr << "OK" <<endl;
					if (coverage < a[ii].second) coverage = a[ii].second;
				}
			}
			DEBUG("seq: s" << s);
			tmp_res.push_back(s);
			tmp_cov.pb(coverage);
//			Sequence* tmpSeq = new Sequence(s);
//			res.pb(make_pair(tmpSeq,coverage));
			color++;
		}
	}
/*	memset(used, 0, sizeof(used));
	vector<ll> lmers[MAXLMERSIZE];
	forn(i, tmp_res.size()) {
		lmers[i].clear;
		char t_seq[200] = tmp_res[i].c_str();
		forn(j, tmp_res[i].length() - l)
			lmers[i].pb(extractMer(t_seq,l, j ));
		sort(lmers[i].begin(), lmers[i].end());
	}*/
	forn(i,  tmp_res.size()) {

		int good = 1;
		for(int j = i + 1; j < tmp_res.size(); j++){
			pair<int, pair<int, int> > comp_res = maxCommonSubstring(tmp_res[i], tmp_res[j]);
			if (comp_res.fi > l) {
				good = 0;
				tmp_res[j] = tmp_res[j].substr(comp_res.se.se, comp_res.fi);
				break;

			}
		}
		if (good) {
			Sequence* tmpSeq = new Sequence(tmp_res[i]);
			res.pb(make_pair(tmpSeq,tmp_cov[i]));

		}
	}
	/*
	{
		forn(i, size) {
			cerr << left[i] << " ";
			cerr << shift_left[i] << " ";

		}
		cerr << endl;
		forn(i, size) {
			cerr << right[i] << " ";
			cerr << shift_right[i] << " ";
		}
		forn(i, res.size()) {
			cerr<<(res[i].first)->str() << endl;
		}
	}
	assert(0);*/
	return res;
}
/*
 * @return (length of intersection, start position in first and second strings)
 *
 */
pair<int, pair<int,int>> maxCommonSubstring(string &s1,string &s2) {
	int l1 = s1.length();
	int l2 = s2.length();
	char table[200][200];
	assert(l1 <200 && l2 <200);
	forn(i, l1)
		forn(j, l2)
			table[i][j] = 0;
	forn(i, l1)
		forn(j, l2)
			if (i> 0 && j > 0 && s1[i] == s2[j])
				table[i][j] = table[i-1][j-1] + 1;
	pair<int, pair<int,int> > res = make_pair(0, make_pair(0,0));
	forn(i, l1)
		forn(j, l2)
			if (table[i][j] > res.fi) {
				res.fi = table[i][j];
				res.se.fi = i - table[i][j];
				res.se.se = j - table[i][j];
			}
}

inline bool checkBoundsForUpper(ll upper) {
	return true;
	if ((upper >= 1<<20) && (upper < upperMax))
		return true;
	else return false;
}

void addPairToTable(myMap& table, ll upper, ll lower) {
	if (table.find(upper) != table.end()) {
		vector<ll>::iterator it = find(table[upper].first.begin(), table[upper].first.end(), lower);
		if (it == table[upper].first.end()) {
			table[upper].first.pb(lower);
			table[upper].second.pb(1);
			++uniqPairs;
		}
		else {
			int index = distance(table[upper].first.begin(), it);
			table[upper].second[index]++;
		}
	} else {
		pair<vector<ll>,vector<int>> tmp;
		tmp.first.clear();
		tmp.first.pb(lower);
		tmp.second.clear();
		tmp.second.pb(1);
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


void constructTable(string inputFile, myMap &table, bool reverse) {
	FILE* inFile = fopen(inputFile.c_str(), "r");
	int count = 0;
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	while (nextReadPair(inFile, upperNuclRead, lowerNuclRead)) {
//		fprintf(stderr, "%s", upperNuclRead);
		if ((strlen(upperNuclRead)<readLength)||(strlen(lowerNuclRead)<readLength)) continue;
		if (reverse) {
			codeRead(upperNuclRead, lowerRead);
			codeRead(lowerNuclRead, upperRead);
		} else {
			codeRead(upperNuclRead, upperRead);
			codeRead(lowerNuclRead, lowerRead);
		}
		processReadPair(table, upperRead, lowerRead);
		if (!(count & (1024*64 - 1)))
			INFO("read number "<<count<<" processed"<<endl);
		count++;
	}
}

void outputTable(string outputFile, myMap &pairedTable) {
	FILE* outFile = fopen(outputFile.c_str(), "w");
	int j = 0;
	for (myMap::iterator iter = pairedTable.begin() ; iter != pairedTable.end(); iter++) {
		pair<ll, pair<vector<ll>, vector<int>>> p = (*iter);
		fprintf(outFile,"%lld %d\n", p.fi, p.se.fi.size());
		forn(i, p.se.fi.size()) {
			fprintf(outFile,"%lld %d ", p.se.fi[i], p.se.se[i]);
		}
		fprintf(outFile, "\n\n");
		if (!(j & (1024*128-1)))
			DEBUG("Pair number" << j << endl);
		j++;
	}
	pairedTable.clear();
	fclose(outFile);
}

void readsToPairs(string inputFile, string outputFile , bool reverse) {

	myMap table;
	INFO("generation of k-l pairs started");
	constructTable(inputFile, table, reverse);
	INFO("generation of k-l pairs finished, dumping to disk");
	outputTable(outputFile, table);
	table.clear();
}
//#define OUTPUT_DECOMPRESSED
int pairsToLmers(string inputFile, string outputFile) {
	FILE* inFile = fopen(inputFile.c_str(), "r");
	FILE* outFile = fopen(outputFile.c_str(), "w");

	cerr<<"pairsToLmers "<<inputFile.c_str()<<"->"<<outputFile.c_str()<<endl;
	int ok = 1;
	ll kmer; int lsize;
	ll lmers[MAXLMERSIZE];
	int covers[MAXLMERSIZE];

	set<ll> lset;
//	set<ll> kset
	int count = 0;
	while (1) {
		count++;
		ok = fscanf(inFile, "%lld %d", &kmer, &lsize);
		if (ok != 2) {
			if (ok > 0) {
				ERROR("error in reads.");
				break;
			}
			else {
				INFO ("Lmers reading finished!!");
				break;
			}
		}
		if (lsize > MAXLMERSIZE) {
			ERROR("TOO MUCH LMERS CORRESPONDING TO ONE k-mer");
			return -2;
		}

		forn(i, lsize) {
			if (fscanf(inFile, "%lld %d", &lmers[i], &covers[i]) != 2) {
				ERROR( "Error in pairsToSequences reading l-mers");
				return -1;
			}
		}
		sort(lmers, lmers + lsize);
		forn(i, lsize) {
			lset.insert(lmers[i]);
		}
	}
	int lsetsize = lset.size();
	fprintf(outFile, "%d\n", lsetsize);
	for(set<ll>::iterator i = lset.begin(); i != lset.end(); i++ ) {
		fprintf(outFile, "%lld ", *i);
	}
	fclose(outFile);
	fclose(inFile);
	return 0;
}

void readLmersSet(string lmerFile, set<long long > & lset)
{
    ll lmersize, tmp;
    FILE *lFile = fopen(lmerFile.c_str(), "r");
    int ok = fscanf(lFile, "%lld", &lmersize);
    if (ok != 1) cerr << "Error in Lmers reading";
    forn(i, lmersize) {
		ok = fscanf(lFile, "%lld", &tmp);
		if (ok != 1) cerr << "Error in Lmers reading";
		lset.insert(tmp);
	}
    fclose(lFile);
}

int pairsToSequences(string inputFile, string lmerFile, string outputFile) {
	FILE* inFile = freopen(inputFile.c_str(), "r", stdin);
    int ok = 1;
    INFO("PairsToSequences started");
    set<ll> lset;
    readLmersSet(lmerFile, lset);
    pair <ll,int> lmers[MAXLMERSIZE];
	ll kmer;
	int lsize;
	FILE* outFile = fopen(outputFile.c_str(), "w");
	int count = 0;

	while (1) {
		count++;
		ok = fscanf(inFile, "%lld %d", &kmer, &lsize);
		if (ok != 2) {
			if (ok > 0) {
				ERROR ("error in reads");
				assert(0);
			}
			else
				INFO ( "Finished!!");
			break;
		}
		if (lsize > MAXLMERSIZE) {
			ERROR ("TO much lmers" ) ;
			return -2;
		}

		forn(i, lsize) {
			if (fscanf(inFile, "%lld %d", &lmers[i].first, &lmers[i].second) != 2) {
				ERROR("Error in pairsToSequences reading l-mers");
				return -1;
			}
		}
		//sort(lmers, lmers + lsize, ComparePairByFirst);
		downSeqs clusters =  clusterize(lmers, lsize, inClusterMaxShift);
//		return 0;
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
			outstring = clusters[i].first->str();
//			assert(outstring.size() >= l);
//			if (outstring.size() >l+2)
//				fprintf(outFile, "%s %d ",outstring.substr(1,outstring.size()-2).c_str(),clusters[i].second);
//			else
				fprintf(outFile, "%s %d ",outstring.c_str(),clusters[i].second);
		}
		fprintf(outFile, "\n");
#ifdef OUTPUT_DECOMPRESSED
		fprintf(decompressed, "\n");
#endif
		//	return 0;
		if (!(count & ((1 << 15) - 1) ))
			DEBUG("k-sequence pairs for k "<< count <<" generated" <<endl);
	}
	INFO("finished");
	fclose(outFile);
	fclose(inFile);
	return 0;
}

#include "common.hpp"
#include "sequence.hpp"
#include "constructHashTable.hpp"
#include "graphio.hpp"
#include "nucl.hpp"
LOGGER("p.constructHashTable");

using namespace std;

typedef vector<pair<Sequence*,int>> downSeqs;


int totalKmers = 0;
int uniqPairs = 0;
int uniqKmers = 0;
const int MAXLMERSIZE = 10000;
#define MAX_COVERAGE 10000000
ll upperMask;
ll lowerMask;

ll upperMax;

map<ll, int> kmers;
map<ll, int> lmers;
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
	char cur_string[MAXLMERSIZE * 2];
	memset(cur_string, 0, sizeof(cur_string));
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
	char cur_string[MAXLMERSIZE * 2];
	vector<int> tmp_cov;
	assert (max_shift <= 20);

	//cerr << "Start clustering"<<endl;
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
				if (used[left[ii]] == color) break;
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
				if (used[right[ii]] == color) break;
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

				}
				coverage += a[ii].second;
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
//	cerr << "Start compressing"<<endl;

	forn(i,  tmp_res.size()) {

		int good = 1;
		/*for(int j = i + 1; j < tmp_res.size(); j++){
			pair<int, pair<int, int> > comp_res = maxCommonSubstring(tmp_res[i], tmp_res[j]);
			if (comp_res.fi > 2*l - 1) {
				good = 0;
				DEBUG(" FOUND intersection length" << comp_res.fi << " on second position " << comp_res.se.se <<" " << tmp_res[i] <<" "<< tmp_res[j]);
				tmp_res[j] = tmp_res[j].substr(comp_res.se.se, comp_res.fi);
				tmp_cov[j] += tmp_cov[i];
				break;

			}
		}*/
		int s_len = tmp_res[i].length();
		forn(ii, s_len)
			cur_string[ii] = tmp_res[i][ii];
		cur_string[s_len] = 0;
		int sum_cov = 0;
		ll tmpKmer = 0;
		forn(ii, s_len - l) {
			tmpKmer = extractMer(cur_string, ii, l);
			forn(j, size) {
				if (a[j].first == tmpKmer) {
					sum_cov +=a[j].second;
					break;
				}
			}
		}
		DEBUG("sum_cov computed for string " << tmp_res[i] << "with s_len" << s_len);
		int ii = 0;
		int cov = 0;
		do {
			cov = 0;
			tmpKmer = extractMer(cur_string, ii, l);
			ii++;
			forn(j, size) {
				if (a[j].first == tmpKmer) {
					cov = a[j].second;
					break;
				}
			}
			assert(cov > 0);
		} while ((ii < s_len - l) && (sum_cov > cov * 2 * (s_len - l + 1)));
		int left_start = ii - 1;
		ii = s_len - l ;
		cov = 0;
		DEBUG("Left fixed");
		do {
			cov = 0;
			tmpKmer = extractMer(cur_string, ii, l);
			ii--;
			forn(j, size) {
				if (a[j].first == tmpKmer) {
					cov = a[j].second;
					break;
				}
			}
			assert(cov > 0);

		} while ((ii >left_start + l) && (sum_cov > cov * 2 * (s_len - l + 1)));
		tmp_res[i] = tmp_res[i].substr(left_start, ii - left_start + 1 + l);
		tmp_cov[i] = 0;
//		int sum_cov = 0;
		tmpKmer = 0;
		s_len = tmp_res[i].length();

		DEBUG("right fixed");
		DEBUG("For string" <<tmp_res[i] << " "<< s_len);

		forn(ii, s_len)
			cur_string[ii] = tmp_res[i][ii];
		cur_string[s_len] = 0;
		forn(ii, s_len - l) {
			tmpKmer = extractMer(cur_string, ii, l);
			forn(j, size) {
				if (a[j].first == tmpKmer) {
					tmp_cov[i] +=a[j].second;
					break;
				}
			}
		}
		DEBUG("cov computed");
	}
	forn(i,  tmp_res.size()) {

		int good = 1;
		for(int j = i + 1; j < tmp_res.size(); j++){
			pair<int, pair<int, int> > comp_res = maxCommonSubstring(tmp_res[i], tmp_res[j]);
//			if (comp_res.fi > l - 1) {
			if (comp_res.fi == min(tmp_res[i].length(), tmp_res[j].length())) {
				good = 0;
				DEBUG(" FOUND intersection length" << comp_res.fi << " on second position " << comp_res.se.se <<" " << tmp_res[i] <<" "<< tmp_res[j]);
				tmp_res[j] = tmp_res[j].substr(comp_res.se.se, comp_res.fi);
				tmp_cov[j] += tmp_cov[i];
				break;

			}
		}
		if (good && (tmp_cov[i] * 1000 /(tmp_res[i].length() - l + 1) > coverage_cutoff)) {
			Sequence* tmpSeq = new Sequence(tmp_res[i]);
			res.pb(make_pair(tmpSeq,tmp_cov[i]));

		}
	}
//	cerr << "Finish clustering"<<endl;

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
	DEBUG("finishd clustering of kmer");
	return res;
}
/*
 * @return (length of intersection, start position in first and second strings)
 *
 */
pair<int, pair<int,int>> maxCommonSubstring(string &s1,string &s2) {
	int l1 = s1.length();
	int l2 = s2.length();
	char table[400][400];
	assert(l1 < 400 && l2 < 400);
	forn(i, l1 +1)
		forn(j, l2 +1)
			table[i][j] = 0;
	forn(i, l1 + 1)
		forn(j, l2 + 1)
			if (i> 0 && j > 0 && s1[i-1] == s2[j-1])
				table[i][j] = table[i-1][j-1] + 1;
	pair<int, pair<int,int> > res = make_pair(0, make_pair(0,0));
	forn(i, l1 + 1)
		forn(j, l2 + 1)
			if (table[i][j] > res.fi) {
				res.fi = table[i][j];
				res.se.fi = i - table[i][j];
				res.se.se = j - table[i][j];
			}
/*	if (res.fi > 25) {
		cerr << res.fi << s1 << " " << s2;
		assert(0);
	}*/
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
		vector<ll>::iterator it = find(table[upper].first.begin(), table[upper].first.end(), lower);
		if (it == table[upper].first.end()) {
			table[upper].first.pb(lower);
			table[upper].second.pb(1);
			++uniqPairs;
			if (!(uniqPairs&((1<<15)-1)))
				INFO("unique pairs "<<uniqPairs);
		}
		else {
			int index = distance(table[upper].first.begin(), it);
			table[upper].second[index]++;
		}
	} else {
		pair<vector<ll>,vector<short>> tmp;
		tmp.first.clear();
		tmp.first.pb(lower);
		tmp.second.clear();
		tmp.second.pb(1);
		table.insert(make_pair(upper, tmp));

		//		cerr<<"inserting"<<table.size();
		++uniqKmers;
		++uniqPairs;
		if (!(uniqPairs&((1<<15)-1)))
			INFO("unique pairs "<<uniqPairs);
		if (!(uniqKmers&((1<<15)-1)))
			INFO("unique Kmers "<<uniqKmers);
	}
}


void processReadPair(myMap& table, char *upperRead, char *lowerRead) {
	int shift = (l - k) / 2;
	int up_len = strlen(upperRead);
	int low_len = strlen(lowerRead);
	int low_shift = readLength - low_len;
	if ((up_len<k)||(low_len<l)) return;
	ll upper = extractMer(upperRead, shift, k);
	ll lower = extractMer(lowerRead, 0, l);
//	cerr <<"\n " <<up_len <<"\n" << low_len;
//	cerr <<(string("\n" )+ upperRead) << (string("\n" )+ lowerRead + "\n");
//	cerr.flush();
//	cerr << "Up_len "<<up_len<<" low_len "<<low_len<<endl;
	ll lowers[MAX_READ_LENGTH+2];
	lowers[0] = lower;
	if (low_len > l)
	forn(j, low_len - l) {
		lower <<= 2;
		lower += codeNucleotide( lowerRead[j + l]);
		if (codeNucleotide( lowerRead[j + l])==-1)
			cerr<<"len "<<low_len<<" pos "<<j<<" in "<<lowerRead;
		lower &= lowerMask;
		lowers[j + 1] = lower;
//		cerr << "lowers "<<j+1<<" "<<lowers[j+1]<<endl;
	}
//	cerr << "lowers_coded"<<endl;
	lower = lowers[0];
	//	fprintf(stderr,"%lld %lld\n", upper, lower);
	int j = 0;
	for (; j + k + shift < up_len; j++) {
		if (checkBoundsForUpper(upper)) {
			for (int jj = max(0, j + low_shift - range_variating); jj < min(low_len - l +1, j + low_shift + range_variating + 1); jj ++)
			{			assert(jj<=low_len-l);
				addPairToTable(table, upper, lowers[jj]);
//				if ((lowers[jj]&3)!=2) assert(0);
			}
			totalKmers++;
		}

		upper <<= 2;
		upper += codeNucleotide(upperRead[j + k + shift]);
		upper &= upperMask;
		if (codeNucleotide( upperRead[j + k + shift])==-1)
			cerr<<"up len "<<up_len<<" pos "<<j<<" in "<<lowerRead;

//		lower <<= 2;
//		lower += lowerRead[j + l];
//		lower &= lowerMask;

		//fprintf(stderr,"%d %d\n", upper, lower);

	}
	if (checkBoundsForUpper(upper)) {
		for (int jj = max(0, j + low_shift - range_variating); jj < min(low_len - l +1, j + low_shift +range_variating + 1); jj ++)
		addPairToTable(table, upper, lowers[jj]);
		totalKmers++;
	}
//	cerr << table.size()<<endl;
}

inline void reverseCompliment(char *upperRead, char* lowerRead){
	int up_len = strlen(upperRead);
	int low_len = strlen(lowerRead);
	char * tmpRead = new char[readLength + 2];
	forn(i, up_len) {
		tmpRead[i] = nucl_complement(upperRead[up_len - 1 - i]);
	}

	forn(i, low_len) {
		upperRead[i] = nucl_complement(lowerRead[low_len - 1 - i]);
	}
	upperRead[low_len] = 0;
	forn(i, up_len) {
		lowerRead[i] = tmpRead[i];
	}
	lowerRead[up_len] = 0;

	delete[] tmpRead;


//
//	forn(i, up_len) {
//		tmpRead[i] = nucl_complement(upperRead[up_len - 1 - i]);
//	}
//	forn(i, up_len) {
//		upperRead[i] = tmpRead[i];
//	}
//
//	forn(i, low_len) {
//		tmpRead[i] = nucl_complement(lowerRead[low_len - 1 - i]);
//	}
//	forn(i, low_len) {
//		lowerRead[i] = tmpRead[i];
//	}
	// cerr << "\n\n" << strlen(upperRead) <<" " << strlen(lowerRead) << "\n";
}
void constructTable(string inputFile, myMap &table, bool reverse) {
	FILE* inFile = fopen(inputFile.c_str(), "r");
	int count = 0;
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char* fictiveRead = new char[readLength + 2];

	forn(i, readLength)
		fictiveRead[i] = 'G';
	fictiveRead[readLength] = 0;

	while (nextReadPair(inFile, upperNuclRead, lowerNuclRead)) {
	//	fprintf(stderr, "%s", upperNuclRead);
		// cerr.flush();
//		if ((strlen(upperNuclRead) < readLength)||(strlen(lowerNuclRead) < readLength)){
//			continue;
//		}

	//	cerr << "?";
	//	cerr.flush();
		forn(tmp, 2) {
			if (fictiveSecondReads) {
				processReadPair(table, upperNuclRead, fictiveRead);
				processReadPair(table, lowerNuclRead, fictiveRead);
			} else {
				processReadPair(table, upperNuclRead, lowerNuclRead);
//				statsReadPair(kmers, lmers, upperNuclRead, lowerNuclRead);
			}
			if (!useRevertedPairs)
				break;
			else {
				// cerr << "?";
		//		cerr.flush();
//				if (fictiveSecondReads) {
//	//				reverseCompliment(upperNuclRead, upperNuclRead);
////					reverseCompliment(lowerNuclRead, lowerNuclRead);
//					reverseCompliment(upperNuclRead, lowerNuclRead);
//						} else
//				{
				reverseCompliment(upperNuclRead, lowerNuclRead);
//				}
			}

		}
		if (!(count & (1024*64 - 1)))
			INFO("read number "<<count<<" processed"<<endl);
		count++;
	}
	cerr << "often kmers, rare lmers:\n";
	for (map<ll, int>::iterator it = kmers.begin(); it !=kmers.end(); it++) {
		if ((it->second > 30) && ((lmers.find(it->first) == lmers.end()) || lmers[it->first] < 5)) {
			cerr << decompress(it->first, k) << " " <<it->second << " ";
			cerr << it->first << " " <<it->second << " ";
			if ((lmers.find(it->first)) == lmers.end()) {
				cerr << 0 << endl;
			} else {
				cerr << lmers[it->first] << endl;
			}
        }
	}
/*
	cerr << "often lmers, rare kmers:\n";
	for (map<ll, int>::iterator it = lmers.begin(); it !=lmers.end(); it++) {
				   if ((it->second > 30) && ((kmers.find(it->first) == kmers.end()) || kmers[it->first] < 5)) {
						cerr << it->first << " " <<it->second << " ";
   					   cerr << decompress(it->first, k) << " " <<it->second << " ";
						   if ((kmers.find(it->first)) == kmers.end()) {
								   cerr << 0 << endl;
						   } else {
								   cerr << kmers[it->first] << endl;
						   }
				   }
		  }
*/
}

void outputTable(string outputFile, myMap &pairedTable) {
	FILE* outFile = fopen(outputFile.c_str(), "w");
	int j = 0;
	for (myMap::iterator iter = pairedTable.begin() ; iter != pairedTable.end(); iter++) {
		pair<ll, pair<vector<ll>, vector<short>>> p = (*iter);
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
#define OUTPUT_DECOMPRESSED
int pairsToLmers(string inputFile, string outputFile) {
	FILE* inFile = fopen(inputFile.c_str(), "r");
	FILE* outFile = fopen(outputFile.c_str(), "w");

	cerr<<"pairsToLmers "<<inputFile.c_str()<<"->"<<outputFile.c_str()<<endl;
	int ok = 1;
	ll kmer; ll lsize;
	ll lmers[MAXLMERSIZE];
	int covers[MAXLMERSIZE];
	ll cover;

	map<ll, int> lset;
//	set<ll> kset
	int count = 0;
	while (1) {
		count++;
		ok = fscanf(inFile, "%lld %lld", &kmer, &lsize);
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
			ERROR("TOO MANY LMERS CORRESPONDING TO ONE k-mer");
			cerr<<"kmer "<<kmer<<" size "<<lsize<<cerr;
			return -2;
		}

		forn(i, lsize) {
			if (fscanf(inFile, "%lld %lld", &lmers[i], &cover) != 2) {
				ERROR( "Error in pairsToLmers reading l-mers");
				return -1;
			}
			if (cover>MAX_COVERAGE) covers[i] = MAX_COVERAGE;
			else (covers[i] = cover);
		}
		forn(i, lsize) {
			if (lset.find(lmers[i])!= lset.end())
				lset.insert(mp(lmers[i], covers[i]));
			else
				lset[lmers[i]] +=covers[i];
		}
	}
	int lsetsize = lset.size();
//	fprintf(outFile, "%d\n", lsetsize);
	for(map<ll, int>::iterator i = lset.begin(); i != lset.end(); i++ ) {
		fprintf(outFile, "%lld %d\n", i->first, i->second);
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
  //  set<ll> lset;
  //  readLmersSet(lmerFile, lset);
    pair <ll,int> lmers[MAXLMERSIZE];
	ll kmer;
	int lsize;
	FILE* outFile = fopen(outputFile.c_str(), "w");
#ifdef OUTPUT_DECOMPRESSED
	FILE* decompressed = fopen((outputFile+".decompr").c_str(), "w");
#endif

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
		downSeqs clusters =  clusterize0704(lmers, lsize, inClusterMaxShift);
//		return 0;
		int clsize = clusters.size();
		string outstring;

//		string s = decompress(kmer, k);
#ifdef OUTPUT_DECOMPRESSED
		fprintf(decompressed, "%s %d\n", decompress(kmer, k).c_str(), lsize);
#endif
		fprintf(outFile, "%lld %d\n", kmer, clsize);
#ifdef OUTPUT_DECOMPRESSED
//		forn(i, lsize) {
//			fprintf(decompressed, "%s ", decompress(lmers[i].first, l).c_str());
//		}
#endif
		forn(i, clsize) {
			outstring = clusters[i].first->str();
//			assert(outstring.size() >= l);
//			if (outstring.size() >l+2)
//				fprintf(outFile, "%s %d ",outstring.substr(1,outstring.size()-2).c_str(),clusters[i].second);
//			else
				fprintf(outFile, "%s %d ",outstring.c_str(),clusters[i].second);
#ifdef OUTPUT_DECOMPRESSED
				fprintf(decompressed, "%s %d ",outstring.c_str(),clusters[i].second);
#endif
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
#ifdef OUTPUT_DECOMPRESSED
	fclose(decompressed);
#endif
	return 0;
}


void statsReadPair(map<ll, int>& kmers,map<ll, int>& lmers, char *upperRead, char *lowerRead) {
		int shift = (l - k) / 2;
        int up_len = strlen(upperRead);
        int low_len = strlen(lowerRead);
        int low_shift = readLength - low_len;
        if ((up_len<k)||(low_len<l)) return;
        ll upper = extractMer(upperRead, shift, k);
        ll lower = extractMer(lowerRead, 0, l);
        if (lmers.find(lower) != lmers.end()) {
               lmers[lower] ++;
		} else {
                lmers.insert(mp(lower, 1));
        }
        if (kmers.find(upper) != kmers.end()) {
        	kmers[upper] ++;
        } else {
            kmers.insert(mp(lower, 1));
         }
        //     cerr <<"\n " <<up_len <<"\n" << low_len;
 //     cerr <<(string("\n" )+ upperRead) << (string("\n" )+ lowerRead + "\n");
 //     cerr.flush();
 //     cerr << "Up_len "<<up_len<<" low_len "<<low_len<<endl;
        ll lowers[MAX_READ_LENGTH+2];
        lowers[0] = lower;
        if (low_len > l)
        forn(j, low_len - l) {
                lower <<= 2;
                lower += codeNucleotide( lowerRead[j + l]);
                if (codeNucleotide( lowerRead[j + l])==-1)
                        cerr<<"len "<<low_len<<" pos "<<j<<" in "<<lowerRead;
                lower &= lowerMask;
                lowers[j + 1] = lower;
                if (lmers.find(lowers[j]) != lmers.end()) {
                       lmers[lowers[j]] ++;
				} else {
                        lmers.insert(mp(lowers[j], 1));
                }
 //             cerr << "lowers "<<j+1<<" "<<lowers[j+1]<<endl;
        }
        lower = upper;
        lowers[0] = lower;
        if (up_len > l)
        forn(j, up_len - l) {
                lower <<= 2;
                lower += codeNucleotide( upperRead[j + l]);
                lower &= lowerMask;
                lowers[j + 1] = lower;
                if (kmers.find(lowers[j]) != kmers.end()) {
                        kmers[lowers[j]] ++;
                } else {
                        kmers.insert(mp(lowers[j], 1));
                }
 //             cerr << "lowers "<<j+1<<" "<<lowers[j+1]<<endl;
        }

}

void constructReversedReadPairs(string inputFile, string outputFile) {
	FILE* inFile = fopen(inputFile.c_str(), "r");
	FILE* outFile = fopen(outputFile.c_str(), "w");
	int count = 0;
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	while (nextReadPair(inFile, upperNuclRead, lowerNuclRead)) {
		reverseCompliment(upperNuclRead, lowerNuclRead);
		if (!(count & (1024*64 - 1)))
			INFO("read number "<<count<<" reverted"<<endl);
		count++;
		fprintf(outFile, "%s %s\n", upperNuclRead, lowerNuclRead);
	}

}

void computeGlobalStats(string kmers_s, string reads, string statsOutputFile) {
	FILE* kmerFile = fopen(kmers_s.c_str(), "r");
	FILE* inFile = fopen(reads.c_str(), "r");
	FILE* outFile = fopen(statsOutputFile.c_str(), "w");
	map<ll, pair<int, int> > kmers;
	ll kmer;
	char * upperRead = new char[readLength +10];
	char * lowerRead = new char[readLength +10];
	int a, b;
	while(fscanf(kmerFile, "%lld %d %d", &kmer, &a, &b) == 3) {
		kmers.insert(mp(kmer, mp(0, 0)));
	}
	DEBUG(kmers.size());
	int count = 0;
	while (nextReadPair(inFile, upperRead, lowerRead)) {
		DEBUG(count);
		count ++;
		forn(i, 2) {
			int shift = (l - k) / 2;
	        int up_len = strlen(upperRead);
	        int low_len = strlen(lowerRead);
	        int low_shift = readLength - low_len;
	        if ((up_len<k)||(low_len<l)) continue;
	        ll upper = extractMer(upperRead, shift, k);
	        ll lower = extractMer(lowerRead, 0, l);
	        if (kmers.find(lower) != kmers.end()) {
	               kmers[lower].second ++;
	        }
	        if (kmers.find(upper) != kmers.end()) {
	        	kmers[upper].first++;
	        }
	        //     cerr <<"\n " <<up_len <<"\n" << low_len;
	 //     cerr <<(string("\n" )+ upperRead) << (string("\n" )+ lowerRead + "\n");
	 //     cerr.flush();
	 //     cerr << "Up_len "<<up_len<<" low_len "<<low_len<<endl;
	        ll lowers[MAX_READ_LENGTH+2];
	        lowers[0] = lower;
	        if (low_len > l)
	        forn(j, low_len - l) {
	                lower <<= 2;
	                lower += codeNucleotide( lowerRead[j + l]);
	                if (codeNucleotide( lowerRead[j + l])==-1)
	                        cerr<<"len "<<low_len<<" pos "<<j<<" in "<<lowerRead;
	                lower &= lowerMask;
	                lowers[j + 1] = lower;
	                if (kmers.find(lowers[j]) != kmers.end()) {
	                       kmers[lowers[j]].second ++;
					}
	 //             cerr << "lowers "<<j+1<<" "<<lowers[j+1]<<endl;
	        }
	        lower = upper;
	        lowers[0] = lower;
	        if (up_len > l)
	        forn(j, up_len - l) {
	                lower <<= 2;
	                lower += codeNucleotide( upperRead[j + l]);
	                lower &= lowerMask;
	                lowers[j + 1] = lower;
	                if (kmers.find(lowers[j]) != kmers.end()) {
	                        kmers[lowers[j]].first ++;
	                }
	 //             cerr << "lowers "<<j+1<<" "<<lowers[j+1]<<endl;
	        }

			reverseCompliment(upperRead, lowerRead);
		}
	}
	for(map<ll, pair<int, int> >::iterator it = kmers.begin(); it != kmers.end(); it ++) {
		fprintf(outFile, "%lld %d %d\n", it->first, it->second.first, it->second.second);
	}
}




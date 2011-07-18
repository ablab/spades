/*
 * defs.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#ifndef DEFS_H
#define DEFS_H


#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>


using namespace std;

template<class T, class U>
ostream & operator << (ostream &out, const map<T,U> & m) {
	typename map<T,U>::const_iterator it;
	for (it = m.begin(); it != m.end(); it++) {
		out << it->first << " -> " << it->second << endl;
	}
	return out;
}

//print vector of vectors
template<class T>
ostream & operator << (ostream & out, const vector<vector<T > > & v) {
	out << v.size() << endl;
	if (v.size() != 0) {
		out << v[0];
		for (size_t i = 1; i < v.size(); i++) out << endl << v[i];
	}
	return out;
}

//read vector of vectors
template<class T>
istream & operator >> (istream & in, vector<vector<T > > & v) {
	int count;
	vector<T> row;
	in >> count;
	v.reserve(count);
	for (int i = 0; i < count; i++) {
		row.clear();
		in >> row;
		v.push_back(row);
	}

	return in;
}

//print vector of T
template<class T>
ostream & operator << (ostream & out, const vector<T> & v) {
	//out << v.size() << '\t';
	if (v.size() != 0) {
		out << v[0];
		for (size_t i = 1; i < v.size(); i++) out << '\t' << v[i];
	}
	return out;
}
//read vector of T
template<class T>
istream & operator >> (istream & in, vector<T> & v) {
	int count;
	T val;
	in >> count;
	v.reserve(count);
	for (size_t i = 0; i < count; i++) {
		in >> val;
		v.push_back(val);
	}
	return in;
}

void flip( bool &val) { val = !val; }

inline int randNumber (int min, int max) { return min + int((((double) rand()) / RAND_MAX ) * (max - min + 1) ); }

int nt2num(char c) {
	if (c=='A') return 0;
	if (c=='C') return 1;
	if (c=='G') return 2;
	if (c=='T') return 3;
	if (c=='a') return 0;
	if (c=='c') return 1;
	if (c=='g') return 2;
	if (c=='t') return 3;
	cerr << "Invalid character " << c << " in nt2num." << endl;
	exit(1);
}
char num2nt(int x) {
	if (x==0) return 'A';
	if (x==1) return 'C';
	if (x==2) return 'G';
	if (x==3) return 'T';
	cerr << "Invalid number " << x << " in num2nt.\n";
	exit(1);
}

void open_file(ifstream & inFile, string filename) {
	inFile.open(filename.c_str());
	if (!inFile) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}

void open_file(ofstream & file, string filename) {
	file.open(filename.c_str());
	if (!file) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}

void open_file_binary(ofstream & file, string filename) {
	file.open(filename.c_str(), ios::out | ios::binary);
	if (!file) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}
string read_genome(string filename) {
	string genome;
	ifstream inFile;
	open_file(inFile, filename);
	char buf[257];

	//read the file into genome
	inFile.get(buf, 256);

	char c;
	while ((c = inFile.rdbuf()->sbumpc()) != EOF) {
		if (c != '\n') {
			genome.push_back(c);
		}
	}
	inFile.close();
	return genome;
}

template<class T>
int argmin(const vector<T> & data) {
	if (data.size() == 0) return -1;
	T maxVal;
	int maxPos = 0;
	maxVal = data[0];
	for (size_t i = 1; i < data.size(); i++) {
		if (data[i] < maxVal) {
			maxVal = data[i];
			maxPos = i;
		}
	}
	return maxPos;
}

template<class T>
int argmin(T data[], int size) {
	if (size == 0) return -1;
	T maxVal;
	int maxPos = 0;
	maxVal = data[0];
	for (int i = 1; i < size; i++) {
		if (data[i] < maxVal) {
			maxVal = data[i];
			maxPos = i;
		}
	}
	return maxPos;
}


template<class T>
int argmax(const vector<T> & data) {
	if (data.size() == 0) return -1;
	T maxVal;
	int maxPos = 0;
	maxVal = data[0];
	for (int i = 1; i < data.size(); i++) {
		if (data[i] > maxVal) {
			maxVal = data[i];
			maxPos = i;
		}
	}
	return maxPos;
}

template<class T>
int argmax(T data[], int size) {
	if (size == 0) return -1;
	T maxVal;
	int maxPos = 0;
	maxVal = data[0];
	for (int i = 1; i < size; i++) {
		if (data[i] > maxVal) {
			maxVal = data[i];
			maxPos = i;
		}
	}
	return maxPos;
}

int nscore(vector<int> & data, int genLen, int nval) {
	sort(data.begin(), data.end());
	int totLen=0;
	for (size_t i = 0; i < data.size(); i++) {
		totLen += data[i];
		if (totLen >= genLen * double(nval / 100.0)) return data[i];
	}
	cerr << "trouble in nscore\n";
	exit(1);
}


bool get_row(istream & inFile, vector<string> & row, char delim = '\t') {  //read a tab delimited line from file
	string line, s;
	row.clear();
	getline(inFile, line);
	if (inFile.eof()) return false;
	istringstream lineStream(line);

	while (getline(lineStream, s, delim)) {
		row.push_back(s);
	}
	return true;
}

bool get_row_whitespace (istream & inFile, vector<string> & row ) {  //read a tab delimited line from file
	string line, s;
	row.clear();
	getline(inFile, line);
	if (inFile.eof()) return false;
	istringstream lineStream(line);

	while (lineStream >> s) {
		row.push_back(s);
	}
	return true;
}

char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
	return 'X';
}

string revcomp (string s) {
	string rc;
	for (int i = s.length() - 1; i >= 0; i--) rc+=revcomp(s[i]);
	return rc;
}

string rcnorm (string &s) {
	string rc;
	for (size_t i = 0; i < s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]); 
		if (s[i] < c) {
			return s;
		} else if (s[i] > c) {
			return revcomp(s);
		}
	}
	return revcomp(s);
}

string rcnorm (string &s, bool & flipped) {
	string rc;
	for (size_t i = 0; i < s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]); 
		if (s[i] < c) {
			flipped = false;
			return s;
		} else if (s[i] > c) {
			flipped = true;
			return revcomp(s);
		}
	}
	return revcomp(s);
}
enum HAMDIST_OPTS  {SAME_STRAND, DIFF_STRAND, ANY_STRAND};

int hamdist(const string & x, const string & y, HAMDIST_OPTS mode, int maxAllowed = -1) {
	if (mode == ANY_STRAND) {
		string rcy = revcomp(y);
		return min(hamdist(x,y, SAME_STRAND), hamdist(x,rcy,SAME_STRAND));
	} else if (mode == DIFF_STRAND) {
		string rcy = revcomp(y);
		return hamdist(x,rcy,SAME_STRAND);
	}
	assert (mode == SAME_STRAND);
	assert (x.length() == y.length());
	int dist = 0;
	if (maxAllowed == -1) maxAllowed = x.length();
	for (size_t i = 0; i < x.length(); i++) {
		if (x[i] != y[i]) {
			dist++;
			if (dist > maxAllowed) return dist;
		}
	}
	return dist;
}

int hamdist(string & x, string & y, HAMDIST_OPTS mode, bool &flipped) {
	assert (mode == ANY_STRAND);
	int dist, dist2;
	dist = hamdist(x,y, SAME_STRAND);
	dist2 = hamdist(x,y,DIFF_STRAND);
	if (dist2 < dist) {
		flipped = true;
		return dist2;
	}
	flipped = false;
	return dist;

}

void maskStrings(vector<string> & strings) {
	size_t minLen = 10000000;
	for (size_t i = 0; i < strings.size(); i++)  {
		if (strings[i].length() < minLen) minLen = strings[i].length();
	}
	for (size_t j = 0; j < minLen; j++) {
		bool same = true;;
		char lastchar = strings[0].at(j);
		for (size_t i = 1; i < strings.size(); i++) {
			if (strings[i].at(j) != lastchar) {
				same = false;
				break;
			}
		}
		if (same) {
			for (size_t i = 0; i < strings.size(); i++) {
				strings[i].at(j) = '-';
			}
		}
	}
}

template<class T>
string make_string(const T & s) {
	ostringstream o;
	o << s;
	return o.str();
}

string add_commas(int num) {
	string s, retval;
	s = make_string(num);
	for (size_t i = 0; i < s.length(); i++) {
		retval.push_back(s.at(i));
		size_t j = s.length() - 1 - i;
		if (((j % 3) == 0) && (j != 0)) {
			retval.push_back(',');
		}
	}
	return retval;
}

#endif




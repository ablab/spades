//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * defs.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

#ifndef DEFS_H
#define DEFS_H


#include <cassert>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <list>
#include <queue>
#include <cstdarg>
#include <algorithm>


using namespace std;




inline int nt2num(char c) {
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
inline char num2nt(int x) {
	if (x==0) return 'A';
	if (x==1) return 'C';
	if (x==2) return 'G';
	if (x==3) return 'T';
	cerr << "Invalid number " << x << " in num2nt.\n";
	exit(1);
}

inline void open_file(ifstream & inFile, string filename) {
	inFile.open(filename.c_str());
	if (!inFile) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}


inline bool get_row(istream & inFile, vector<string> & row, char delim = '\t') {  //read a tab delimited line from file
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

inline bool get_row_whitespace (istream & inFile, vector<string> & row ) {  //read a tab delimited line from file
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

inline char revcomp (char s) {
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

inline string revcomp (string s) {
	string rc;
	for (int i = s.length() - 1; i >= 0; i--) rc+=revcomp(s[i]);
	return rc;
}

inline string rcnorm (string &s) {
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

inline string rcnorm (string &s, bool & flipped) {
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


#endif




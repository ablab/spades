/*
 * splinter.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */

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
#include<cassert>
#include "defs.hpp"

using namespace std;
ostringstream oMsg;

string encode3toabyte (const string & s)  {
	string retval;
	char c = 48;
	int weight = 16;
	size_t i;
	for (i = 0; i < s.length(); i += 1) {
		if (i % 3 == 0) {
			c= 48;
			weight = 16;
		}
		c += weight * nt2num(s[i]);
		weight /= 4;
		if (i % 3 == 2) retval += c;
	}
	if (i % 3 != 0) retval += c;
	return retval;
}


int main(int argc, char * argv[]) {
	int index = atoi(argv[1]);  
	int modval = atoi(argv[2]);  
	vector<string> row;
	int counter = 0;
	while (get_row_whitespace(cin, row)) {
		string s = row[0];
		string src = revcomp(s);
		string sub;
		for (size_t i = index; i < s.length(); i += modval) {
			sub += src.at(i);
		}
		//cout << sub << "\t" << counter << endl;
		cout << encode3toabyte(sub) << "\t" << counter << endl;
		sub = "";
		for (size_t i = index; i < s.length(); i += modval) {
			sub += s.at(i);
		}
		//cout << sub << "\t" << counter << endl;
		cout << encode3toabyte(sub) << "\t" << counter << endl;
		counter++;
	}
}

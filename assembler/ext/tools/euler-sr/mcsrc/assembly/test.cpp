/***************************************************************************
 * Title:          test.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "IntegralTupleStatic.h"
using namespace std;
int IntegralTuple::tupleSize = 4;
class Bfc {
public:
	ssize_t a: 10;
	ssize_t b: 10;
	ssize_t c: 12;
};

int main() {
	
<<<<<<< .mine
	Bfc bfc, bfc2;
	bfc.a = 7;
	bfc.b = 28;
	bfc.c = 56;

	ofstream out;
	out.open("test.out", std::ios::binary);
	out.write((char*) &bfc, sizeof(Bfc));
	out.close();

	ifstream in;
	in.open("test.out", std::ios::binary);
	in.read((char*) &bfc2, sizeof(Bfc));
	
	cout << bfc2.a << " " << bfc2.b << " " << bfc2.c << endl;
	
=======
	IntegralTuple tuple;
	string str = "actg";
	tuple.StringToTuple((unsigned char*) str.c_str());
	cout.setf(std::ios::hex);
	cout << "tuple " << tuple.tuple << endl;
	string strrep;
	tuple.ToString(strrep);
	cout << "string: " << strrep << endl;
>>>>>>> .r95
	return 0;
}

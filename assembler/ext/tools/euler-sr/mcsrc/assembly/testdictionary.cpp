/***************************************************************************
 * Title:          testdictionary.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntegralTupleStatic.h"

using namespace std;


int main(int argc, char* argv[]) {
	cout << "argc: " << argc << endl;
	if (argc < 4) {
		cout << "usage: testdict listName tupleSize mode (mode=0 no dictionary, mode=1 dictionary)" << endl;
		exit(0);
	}
	std::string listName = argv[1];
	int tupleSize = atoi(argv[2]);
	ssize_t mode = atoi(argv[3]);
	CountedIntegralTuple *list;
	ssize_t listLength;
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);
	cout << "reading " <<  listName << " " << mode << endl;
	DictionaryTupleList<CountedIntegralTuple, 9> dictTupleList;

	//	dictTupleList.InitFromFile(listName);
	ReadBinaryTupleList(listName, &list, listLength);

	cout << "searching for " << listLength << " tuples." << endl;
	ssize_t i;
	ssize_t nFound = 0;
	if (mode == 0) {
		for (i = 0 ; i < listLength; i++ ){
			if (LookupBinaryTuple(list, listLength, list[i]) != -1) {
				++nFound;
			}
		}
	}
	if (mode == 1) {
		dictTupleList.list = list;
		dictTupleList.listLength = listLength;
		dictTupleList.IndexList();
		for (i = 0; i < dictTupleList.listLength; i++ ){
			if (dictTupleList.DictLookupBinaryTuple(dictTupleList.list[i]) != -1) {
				++nFound;
			}
		}
	}
	cout << "found: " << nFound << " tuples " << endl;
}

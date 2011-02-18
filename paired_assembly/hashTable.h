#ifndef HASHTABLE_H_
#define HASHTABLE_H_

struct hashTable;// index: kmer, value: list of second parts

*hashTable createHashTable();

//0 if success
//1 otherwise
int put(hashTable table, kmerType upper, lmerType lower);



#endif /*HASHTABLE_H_*/


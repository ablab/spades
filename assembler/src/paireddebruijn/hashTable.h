#ifndef HASHTABLE_H_
#define HASHTABLE_H_
#include "pairedGraph.h"

class ListItem {
	ListItem *next;
};

class Record {
	Sequence *lower;
	ListItem *next;
};
class HashTable {
public:
	HashTable();
	ListItem *get(Kmer kmer);
	ListItem *put(Kmer kmer, Sequence *pair, ListItem *next);
};

#endif /*HASHTABLE_H_*/


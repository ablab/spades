#ifndef READSITERATOR_H_
#define READSITERATOR_H_

*Iterator createIterator(FILE *input);

void resetIterator();

ReadPair nextRead(Iterator *iterator);

#endif /*READSITERATOR_H_*/

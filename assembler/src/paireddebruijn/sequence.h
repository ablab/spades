/*
 * sequence.h
 *
 *  Created on: 22.02.2011
 *      Author: student
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

namespace paired_assembler {
typedef long long Kmer;

class Sequence {
	char *_nucleotides;
	short _length;
public:
	Sequence(char *nucleotides, short length);
	char operator[](const int &index);
};

}
#endif /* SEQUENCE_H_ */

/*
 * ifastqstreamTest.cpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "ifaststream.hpp"

void TestIFastaStreamNoFile() {
	ifaststream ifs("./no-file");
	ASSERT(!ifs.is_open());
}

void TestIFastaStreamSingleRead() {
	ifaststream ifs("./data/test/s_6_1.fastq.gz");
	string name, seq, qual;
	ifs >> name >> seq >> qual;
	ASSERT_EQUAL("EAS20_8_6_1_2_768/1", name);
	ASSERT_EQUAL("CAGCACAGAGGATATCGCTGTTACANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", seq);
	ASSERT_EQUAL("HGHIHHHGHECHHHHHHHGGHHHHH###########################################################################", qual);
	ifs >> name >> seq >> qual;
	ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", name);
	ASSERT_EQUAL("CTTGGTGCGGAACTGAAAAGTGGTANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", seq);
	ASSERT_EQUAL("GGGGCGGGGEGGGGGBGAGF:CCCC###########################################################################", qual);
}

void TestIFastaStreamFull() {
	ifaststream ifs("./data/test/s_6_1.fastq.gz");
	string name, seq, qual;
	while (!ifs.eof()) {
		ifs >> name >> seq >> qual;
	}
	ifs.close();
	ASSERT_EQUAL("TEST/1", name);
	ASSERT_EQUAL("CATACGGGTTTCCGCCAGTNTTTCCATGCCGCGATGGACGTAGAACAGACGGTAGTCGGCGTCGATAATGTTTTCGCCATCGACGCACAGACGGAAGTGG", seq);
	ASSERT_EQUAL("HHHGHHGIHIHHEHHHHHGHHHHHHHHGHEHHHHDHAHHHA?HFHEFHEHHHGHGGHGG@B2BEBEF=HHGEEA:C?CCD?B@EF/4=2<4188.?BA5=", qual);
}

cute::suite IFastaStreamSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIFastaStreamNoFile));
	s.push_back(CUTE(TestIFastaStreamSingleRead));
	s.push_back(CUTE(TestIFastaStreamFull));
	return s;
}

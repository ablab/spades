/*
 * ifastqstreamTest.cpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "ireadstream.hpp"

void TestIReadStreamNoFile() {
	ireadstream ifs("./no-file");
	ASSERT(!ifs.is_open());
}

void TestIReadStreamSingleRead() {
	ireadstream ifs("./test/data/s_6_1.fastq.gz");
	ASSERT(ifs.is_open());
	Read r;
	ifs >> r;
	ASSERT_EQUAL("EAS20_8_6_1_2_768/1", r.getName());
	ASSERT_EQUAL("CAGCACAGAGGATATCGCTGTTACANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", r.getSequenceString());
	//ASSERT_EQUAL("HGHIHHHGHECHHHHHHHGGHHHHH###########################################################################", r.getQuality());
	ifs >> r;
	ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", r.getName());
	ASSERT_EQUAL("CTTGGTGCGGAACTGAAAAGTGGTANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", r.getSequenceString());
	//ASSERT_EQUAL("GGGGCGGGGEGGGGGBGAGF:CCCC###########################################################################", r.getQuality());
}

void TestIReadStreamFull() {
	ireadstream ifs("./test/data/s_6_1.fastq.gz");
	ASSERT(ifs.is_open());
	Read r;
	while (!ifs.eof()) {
		ifs >> r;
	}
	ifs.close();
	ASSERT_EQUAL("TEST/1", r.getName());
	ASSERT_EQUAL("CATACGGGTTTCCGCCAGTNTTTCCATGCCGCGATGGACGTAGAACAGACGGTAGTCGGCGTCGATAATGTTTTCGCCATCGACGCACAGACGGAAGTGG", r.getSequenceString());
	//ASSERT_EQUAL("HHHGHHGIHIHHEHHHHHGHHHHHHHHGHEHHHHDHAHHHA?HFHEFHEHHHGHGGHGG@B2BEBEF=HHGEEA:C?CCD?B@EF/4=2<4188.?BA5=", r.getQuality());
}

cute::suite IReadStreamSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIReadStreamNoFile));
	s.push_back(CUTE(TestIReadStreamSingleRead));
	s.push_back(CUTE(TestIReadStreamFull));
	return s;
}

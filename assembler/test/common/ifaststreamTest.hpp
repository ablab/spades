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
	ifaststream ifs("./data/short/s_6_1.fastq.gz");
	string name, seq, qual;
	ifs >> name >> seq >> qual;
	ASSERT_EQUAL("EAS20_8_6_1_2_768/1", name);
	ASSERT_EQUAL("CAGCACAGAGGATATCGCTGTTACANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", seq);
	ASSERT_EQUAL("HGHIHHHGHECHHHHHHHGGHHHHH###########################################################################", qual);
}

cute::suite IFastaStreamSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIFastaStreamNoFile));
	s.push_back(CUTE(TestIFastaStreamSingleRead));
	return s;
}

/*
 * parserTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "ireadstream.hpp"

void TestIReadStream() {
	ireadstream<100,2,char> irs("./data/short/s_6_1.fastq.gz", "./data/short/s_6_2.fastq.gz");
	strobe_read<100,2,char> sr;
	mate_read<100,char>::type mr;
	irs >> sr >> mr;
	ASSERT_EQUAL("CATTATTAGGGATGATTGTGACCCGCGTCAGACCAATCAAATTCGCCAGCGTTTCCACGGGTTTTAGATGACCATAGTGCACCGGATCAAAGGTGCCGCC", sr.get(0).str());
	ASSERT_EQUAL("ATTACGGTCAGTCAGTGTGGGCAGAGCTGGAAGGGTTATCTCTTCTGTTGTGCCATAAACCCCTGGCGGACGTATTTATCGACGGTTGATATGTAATCTT", sr.get(1).str());
	ASSERT_EQUAL("GCAGTAAAGCTATCATGGCAGAATCATTTGCAACTGGTTCCGATCATCAGGTTGTAAACGAGCTCAACGGGGAAAGACTGAGAGAACCAAACGACGTTTT", mr.get(0).str());
	ASSERT_EQUAL("GGCCAACGTTAATTTTGTTACCGACTAAAGTAGAAACTATTTCTTTTAGATGGTCGCATTTATATTTTGCATCGTCCACTTGAAAATCATATCTTATTGC", mr.get(1).str());
}

void TestIReadStreamFull() {
	ireadstream<100,2,char> irs("./data/short/s_6_1.fastq.gz", "./data/short/s_6_2.fastq.gz");
	mate_read<100,char>::type mr;
	while (!irs.eof()) {
		irs >> mr;
	}
	irs.close();
	ASSERT_EQUAL("CATACGGGTTTCCGCCAGTTTTTCCATGCCGCGATGGACGTAGAACAGACGGTAGTCGGCGTCGATAATGTTTTCGCCATCGACGCACAGACGGAAGTGG", mr.get(0).str());
	ASSERT_EQUAL("CGTCCGGCACCGACCACCGATGCTGAAACCTACGAGTTCATCAACGAACTGGGCGACAAGAAAAACAACGTCGTGCCGATTGGTCCGCTGCACGTCACTT", mr.get(1).str());

}

cute::suite IReadStreamSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIReadStream));
	s.push_back(CUTE(TestIReadStreamFull));
	return s;
}

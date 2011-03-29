/*
 * parserTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "ireadstream.hpp"
#include "quality_read_stream.hpp"

static string files[] = {"./test/data/s_6_1.fastq.gz", "./test/data/s_6_2.fastq.gz"};

void TestIReadStream() {
	ireadstream<100,2,short> irs(files);
	strobe_read<100,2,short> sr;
	mate_read<100,short>::type mr;
	irs >> sr >> mr;
	ASSERT_EQUAL("CATTATTAGGGATGATTGTGACCCGCGTCAGACCAATCAAATTCGCCAGCGTTTCCACGGGTTTTAGATGACCATAGTGCACCGGATCAAAGGTGCCGCC", sr[0].str());
	ASSERT_EQUAL("ATTACGGTCAGTCAGTGTGGGCAGAGCTGGAAGGGTTATCTCTTCTGTTGTGCCATAAACCCCTGGCGGACGTATTTATCGACGGTTGATATGTAATCTT", sr[1].str());
	ASSERT_EQUAL("GCAGTAAAGCTATCATGGCAGAATCATTTGCAACTGGTTCCGATCATCAGGTTGTAAACGAGCTCAACGGGGAAAGACTGAGAGAACCAAACGACGTTTT", mr[0].str());
	ASSERT_EQUAL("GGCCAACGTTAATTTTGTTACCGACTAAAGTAGAAACTATTTCTTTTAGATGGTCGCATTTATATTTTGCATCGTCCACTTGAAAATCATATCTTATTGC", mr[1].str());
}

void TestIReadStreamFull() {
	ireadstream<100,2,short> irs(files);
	mate_read<100,short>::type mr;
	while (!irs.eof()) {
		irs >> mr;
	}
	irs.close();
	ASSERT_EQUAL("CATACGGGTTTCCGCCAGTTTTTCCATGCCGCGATGGACGTAGAACAGACGGTAGTCGGCGTCGATAATGTTTTCGCCATCGACGCACAGACGGAAGTGG", mr[0].str());
	ASSERT_EQUAL("CGTCCGGCACCGACCACCGATGCTGAAACCTACGAGTTCATCAACGAACTGGGCGACAAGAAAAACAACGTCGTGCCGATTGGTCCGCTGCACGTCACTT", mr[1].str());
}

void TestIReadStreamReset() {
	ireadstream<100,2,short> irs(files);
	strobe_read<100,2,short> sr;
	mate_read<100,short>::type mr;
	while (!irs.eof()) {
		irs >> mr;
	}
	irs.reset();
	irs >> sr >> mr;
	ASSERT_EQUAL("CATTATTAGGGATGATTGTGACCCGCGTCAGACCAATCAAATTCGCCAGCGTTTCCACGGGTTTTAGATGACCATAGTGCACCGGATCAAAGGTGCCGCC", sr[0].str());
	ASSERT_EQUAL("ATTACGGTCAGTCAGTGTGGGCAGAGCTGGAAGGGTTATCTCTTCTGTTGTGCCATAAACCCCTGGCGGACGTATTTATCGACGGTTGATATGTAATCTT", sr[1].str());
	ASSERT_EQUAL("GCAGTAAAGCTATCATGGCAGAATCATTTGCAACTGGTTCCGATCATCAGGTTGTAAACGAGCTCAACGGGGAAAGACTGAGAGAACCAAACGACGTTTT", mr[0].str());
	ASSERT_EQUAL("GGCCAACGTTAATTTTGTTACCGACTAAAGTAGAAACTATTTCTTTTAGATGGTCGCATTTATATTTTGCATCGTCCACTTGAAAATCATATCTTATTGC", mr[1].str());
}

//void TestQuality1() {
//	QualityReadStream qrs("./test/data/s_6_1.fastq.gz");
//	while (!qrs.eof()) {
//		pair<Sequence, vector<int> > pair = qrs.Next();
//		Sequence s = pair.first;
//		vector<int> q = pair.second;
//		cout << s.str() << endl;
//		for (size_t i = 0; i < q.size(); ++i) {
//			cout << q[i] << " ";
//		}
//		cout << endl;
//	}
//}

cute::suite IReadStreamSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIReadStream));
	s.push_back(CUTE(TestIReadStreamFull));
	s.push_back(CUTE(TestIReadStreamReset));
//	s.push_back(CUTE(TestQuality1));
	return s;
}

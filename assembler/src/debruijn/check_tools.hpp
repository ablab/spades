//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * check_tools.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef CHECK_TOOLS_HPP_
#define CHECK_TOOLS_HPP_

#include "omni/paired_info.hpp"
#include "new_debruijn.hpp"
#include "graph_pack.hpp"

namespace debruijn_graph {


/*void FindZeros(PairedInfoIndex<Graph>& etalon_paired_index) {
	for (auto it = etalon_paired_index.begin(); it != etalon_paired_index.end(); ++it) {
		const vector<PairInfo<EdgeId>> infos = *it;
		for (auto it2 = infos.begin(); it2!=infos.end(); ++it2) {
			PairInfo<EdgeId> info = *it2;
			if (info.first == info.second && info.d == 0.) {
				cout << "FOUND ZEROS!!!" << endl;
				return;
			}
		}
	}
}*/

//void ToSet(PairedInfoIndex<Graph>& paired_index, set<PairInfo<EdgeId>>& as_set) {
////	static size_t count = 0;
//	for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
//		vector<PairInfo<EdgeId>> infos = *it;
//		for (auto it2 = infos.begin(); it2!=infos.end(); it2++) {
////			count++;
//			as_set.insert(*it2);
////			cout << "HERE1 " << *it2 << endl;
//		}
//	}
////	cout << "Count was " << count << endl;
//}

void CheckPairInfo(const vector<PairInfo<EdgeId>>& infos1, const vector<PairInfo<EdgeId>>& infos2) {
  using omnigraph::operator<<;

	for (auto it = infos1.begin(); it != infos1.end(); ++it) {
		bool found = false;
		for (auto it2 = infos2.begin(); it2 != infos2.end(); ++it2) {
			if (*it == *it2)
				found = true;
		}
		if (!found && !((*it).first == (*it).second && (*it).d == 0.)) {
			cerr << "Didn't find " << *it << " in " << infos2 <<  " initial:" << infos1 << endl;
		}
	}

}

void CheckInfoEquality(PairedInfoIndex<Graph>& paired_index1, PairedInfoIndex<Graph>& paired_index2) {
	for (auto it = paired_index1.begin(); it != paired_index1.end(); ++it) {
		const vector<PairInfo<EdgeId>> infos1 = *it;
		EdgeId first = infos1.front().first;
		EdgeId second = infos1.front().second;
		const vector<PairInfo<EdgeId>> infos2 = paired_index2.GetEdgePairInfo(first, second);
		CheckPairInfo(infos1, infos2);
		CheckPairInfo(infos2, infos1);
	}
//	set<PairInfo<EdgeId>> set1;
//	set<PairInfo<EdgeId>> set2;
//	ToSet(paired_index1, set1);
//	ToSet(paired_index2, set2);
//	if (set1.size() != set2.size()) {
//		cerr << "HELP_1!!!" << endl;
//	}
//	for (auto it = set1.begin(); it != set1.end(); ++it) {
//		if (set2.count(*it) == 0) {
//			cerr << "HELP_2!!!" << endl;
//		}
//	}
}

bool CheckContains(const runtime_k::RtSeq& pattern, const Sequence& s, size_t k) {
	if (s.size() < k)
		return false;

	runtime_k::RtSeq kmer = s.start<runtime_k::RtSeq::max_size>(k);
	kmer >>= 0;
	for (size_t i = k - 1; i < s.size(); ++i) {
		kmer <<= s[i];
		if (pattern == kmer)
			return true;
	}
	return false;
}

bool CheckContainsSubKmer(const Sequence& pattern, const Sequence& s, size_t k) {
	if (pattern.size() < k || s.size() < k)
		return false;

	runtime_k::RtSeq kmer = pattern.start<runtime_k::RtSeq::max_size>(k);
	kmer >>= 0;
	for (size_t i = k - 1; i < pattern.size(); ++i) {
		kmer <<= pattern[i];
		if (CheckContains(kmer, s, k)){
			//cout << "Kmer " << kmer << endl;
            return true;
        }
	}
	return false;
}

size_t ThreadedPairedReadCount(const Sequence& s1, const Sequence& s2, io::IReader<io::PairedRead>& stream, size_t k) {
	size_t count = 0;
	io::PairedRead paired_read;
	while (!stream.eof()) {
		stream >> paired_read;
		Sequence read_s1 = paired_read.first().sequence();
		Sequence read_s2 = paired_read.second().sequence();
		if ((CheckContainsSubKmer(read_s1, s1, k) && CheckContainsSubKmer(read_s2, s2, k))
           || (CheckContainsSubKmer(read_s1, s2, k) && CheckContainsSubKmer(read_s2, s1, k))
            ) {
			count++;
            //cout << "Read first " << read_s1 << endl 
                 //<< "Seq  first " << s1 << endl 
                 //<< "Read  second " << read_s2 << endl
                 //<< "Seq   second " << s2 << endl;
		}
	}
	return count;
}

size_t ThreadedPairedReadCount(const conj_graph_pack& gp, int e1, int e2, io::IReader<io::PairedRead>& stream, size_t k) {
	return ThreadedPairedReadCount(gp.g.EdgeNucls(gp.int_ids.ReturnEdgeId(e1)), gp.g.EdgeNucls(gp.int_ids.ReturnEdgeId(e2)), stream, k);
}

double TotalPositiveWeight(const conj_graph_pack& gp, PairedInfoIndex<Graph> paired_index, int e1, int e2) {
	vector<PairInfo<EdgeId>> infos = paired_index.GetEdgePairInfo(gp.int_ids.ReturnEdgeId(e1), gp.int_ids.ReturnEdgeId(e2));
	double s = 0.;
	for (auto it = infos.begin(); it != infos.end(); ++it) {
		double weight = it->weight;
        TRACE("Weight " << weight);
		//if (math::gr(weight, 0.)) {
			s += weight;
		//}
	}
	return s;
}

}

#endif /* CHECK_TOOLS_HPP_ */

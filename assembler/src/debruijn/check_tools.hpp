/*
 * check_tools.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef CHECK_TOOLS_HPP_
#define CHECK_TOOLS_HPP_

#include "paired_info.hpp"
#include "new_debruijn.hpp"

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

void ToSet(PairedInfoIndex<Graph>& paired_index, set<PairInfo<EdgeId>>& as_set) {
	for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
		vector<PairInfo<EdgeId>> infos = *it;
		for (auto it2 = infos.begin(); it2!=infos.end(); it2++) {
			as_set.insert(*it2);
		}
	}
}

void CheckPairInfo(const vector<PairInfo<EdgeId>>& infos1, const vector<PairInfo<EdgeId>>& infos2) {

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
//	//todo remove assert
//	if (set1.size() != set2.size()) {
//		cerr << "HELP_1!!!" << endl;
//	}
//	for (auto it = set1.begin(); it != set1.end(); ++it) {
//		if (set2.count(*it) == 0) {
//			cerr << "HELP_2!!!" << endl;
//		}
//	}
}

}

#endif /* CHECK_TOOLS_HPP_ */

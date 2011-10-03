/*
 * paired_info_checker.cpp
 *
 *  Created on: Sep 26, 2011
 *      Author: andrey
 */

#include "../lc_common.hpp"
#include "../lc_io.hpp"

class PairedInfoChecker {
private:
	Graph& g_;

public:
	PairedInfoChecker(Graph& g) : g_(g) {

	}

	bool IsSymmetric(PairedInfoIndex<Graph>& index) {
		bool result = true;
		for (auto iter = index.begin(); iter != index.end(); ++iter) {
			auto pi = *iter;
			if (pi.size() == 0) {
				continue;
			}
			EdgeId e1 = pi.back().first;
			EdgeId e2 = pi.back().second;

			auto sym_pi = index.GetEdgePairInfo(e2, e1);

			for (auto i1 = pi.begin(); i1 != pi.end(); ++i1) {
				for (auto i2 = sym_pi.begin(); i2 != sym_pi.end(); ++i2) {
					if (math::eq(i1->d, - i2->d) && !math::eq(i1->weight, i2->weight)) {
						INFO("No symmetric found ");
						result = false;
					}
				}
			}

		}
		return result;
	}

	bool IsConjugateSymmetric(PairedInfoIndex<Graph>& index) {
		bool result = true;
		for (auto iter = index.begin(); iter != index.end(); ++iter) {
			auto pi = *iter;
			if (pi.size() == 0) {
				continue;
			}
			EdgeId e1 = pi.back().first;
			EdgeId e2 = pi.back().second;

			auto conj_pi = index.GetEdgePairInfo(g_.conjugate(e1), g_.conjugate(e2));

			for (auto i1 = pi.begin(); i1 != pi.end(); ++i1) {
				for (auto i2 = conj_pi.begin(); i2 != conj_pi.end(); ++i2) {
					double new_d = i1->d - g_.length(e1) + g_.length(e2);
					if (math::eq(i1->d, - new_d) && !math::eq(i1->weight, i2->weight)) {
						INFO("No conjugate found ");
						result = false;
					}
				}
			}

		}
		return result;
	}

	bool AreEqual(PairedInfoIndex<Graph>& index1, PairedInfoIndex<Graph>& index2) {
		bool result = true;
		for (auto iter = index1.begin(); iter != index1.end(); ++iter) {
			auto pi = *iter;
			if (pi.size() == 0) {
				continue;
			}
			EdgeId e1 = pi.back().first;
			EdgeId e2 = pi.back().second;

			auto pi2 = index2.GetEdgePairInfo(e1, e2);

			for (auto i1 = pi.begin(); i1 != pi.end(); ++i1) {
				for (auto i2 = pi2.begin(); i2 != pi2.end(); ++i2) {
					if (math::eq(i1->d, i2->d) && !math::eq(i1->weight, i2->weight)) {
						INFO("Unequal weights");
						result = false;
					}
				}
			}

		}
		return result;

		return false;
	}
};


int main() {
	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	IdTrackHandler<Graph> intIds(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	KmerMapper<K+1, Graph> mapper(g);
	Sequence sequence("");

	long_contigs::LoadFromFile(lc_cfg::get().ds.graph_file, &g, &intIds, sequence, &mapper);
	PairedInfoChecker checker(g);

	DataScanner<Graph> dataScanner(g, intIds);
	dataScanner.loadPaired(lc_cfg::get().u.file1, pairedIndex);

	switch (lc_cfg::get().u.mode) {
	case 1: {
		INFO("Checking " << lc_cfg::get().u.file1);
		INFO("Symmetric: " << checker.IsSymmetric(pairedIndex));
		INFO("Conjugate symmetric: " << checker.IsConjugateSymmetric(pairedIndex));
		break;
	}
	case 2: {
		PairedInfoIndex<Graph> pairedIndex2(g, 0);
		dataScanner.loadPaired(lc_cfg::get().u.file2, pairedIndex2);

		INFO("Checking " << lc_cfg::get().u.file1 << " and " << lc_cfg::get().u.file2);
		INFO("1 is subset of 2 " << checker.AreEqual(pairedIndex, pairedIndex2));
		INFO("2 is subset of 1 " << checker.AreEqual(pairedIndex2, pairedIndex));
		break;
	}
	default: {
		INFO("Unknown mode");
	}
	}

	return 0;
}


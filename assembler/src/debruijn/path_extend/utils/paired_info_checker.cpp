//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * paired_info_checker.cpp
 *
 *  Created on: Sep 26, 2011
 *      Author: andrey
 */

#include "../lc_common.hpp"
#include "../lc_io.hpp"

using namespace debruijn_graph;

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
	}

	void AggregatePairedInfo(PairedInfoIndex<Graph>& clustered, PairedInfoIndex<Graph>& advanced,
			size_t insert_size, size_t read_length,
			PairedInfoIndex<Graph>* result) {

		PairedInfoWeightNormalizer<Graph> normalizer(g_, insert_size, read_length, K);

		for (auto iter = clustered.begin(); iter != clustered.end(); ++iter) {
			auto pi = *iter;
			if (pi.size() == 0) {
				continue;
			}

			EdgeId e1 = pi.back().first;
			EdgeId e2 = pi.back().second;

			auto pi2 = advanced.GetEdgePairInfo(e1, e2);

			for (auto i1 = pi.begin(); i1 != pi.end(); ++i1) {

				auto norm_pi = normalizer.NormalizeWeight(*i1);

				for (auto i2 = pi2.begin(); i2 != pi2.end(); ++i2) {
					if (math::ge(i1->d, i2->d - lc_cfg::get().u.dev) && math::le(i1->d, i2->d + lc_cfg::get().u.dev) && math::gr(i2->weight, 0.0)) {
						norm_pi.weight *= lc_cfg::get().es.advanced_coeff;
					}
				}

				result->AddPairInfo(norm_pi, false);
			}

		}

	}

};


int main() {
	cfg::create_instance(cfg_filename);
	lc_cfg::create_instance(long_contigs::lc_cfg_filename);

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	KmerMapper<K+1, Graph> mapper(g);
	Sequence sequence("");

	long_contigs::LoadFromFile(lc_cfg::get().ds.graph_file, &g, sequence, &mapper);
	PairedInfoChecker checker(g);

	DataScanner<Graph> dataScanner(g);

	switch (lc_cfg::get().u.mode) {
	case 1: {
		INFO("Checking " << lc_cfg::get().u.file1);
		dataScanner.loadPaired(lc_cfg::get().u.file1, pairedIndex);
		INFO("Symmetric: " << checker.IsSymmetric(pairedIndex));
		INFO("Conjugate symmetric: " << checker.IsConjugateSymmetric(pairedIndex));
		break;
	}
	case 2: {
		PairedInfoIndex<Graph> pairedIndex2(g, 0);
		dataScanner.loadPaired(lc_cfg::get().u.file1, pairedIndex);
		dataScanner.loadPaired(lc_cfg::get().u.file2, pairedIndex2);

		INFO("Checking " << lc_cfg::get().u.file1 << " and " << lc_cfg::get().u.file2);
		INFO("1 is subset of 2 " << checker.AreEqual(pairedIndex, pairedIndex2));
		INFO("2 is subset of 1 " << checker.AreEqual(pairedIndex2, pairedIndex));
		break;
	}
	case 3: {
		INFO("Aggregating paired info");

		PairedInfoIndex<Graph> cl(g, 0);
		PairedInfoIndex<Graph> ad(g, 0);
		PairedInfoIndex<Graph> res(g, 0);

		dataScanner.loadPaired(lc_cfg::get().u.clustered, cl);
		dataScanner.loadPaired(lc_cfg::get().u.advanced, ad);

		checker.AggregatePairedInfo(cl, ad,
				lc_cfg::get().u.insert_size, lc_cfg::get().u.read_size,
				&res);

		DataPrinter<Graph> dataPrinter(g);
		dataPrinter.savePaired( "./" + lc_cfg::get().paired_info_file_prefix + "IS" + ToString(lc_cfg::get().u.insert_size) + "_RS" + ToString(lc_cfg::get().u.read_size)
				+ "_agregate_" + ToString(lc_cfg::get().es.advanced_coeff), res);

		INFO("Done");
		break;

	}
	default: {
		INFO("Unknown mode");
	}
	}

	return 0;
}


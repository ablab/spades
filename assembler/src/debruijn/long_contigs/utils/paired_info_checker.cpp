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
					if (i1->d == - i2->d && i1->weight != i2->weight) {
						INFO("No symmetric found: ");
						result = false;
					}
				}
			}

		}
		return result;
	}



	bool AreEqual(PairedInfoIndex<Graph>& index1, PairedInfoIndex<Graph>& index2) {
		return false;
	}
};


int main(int argc, char** argv) {
	if (argc != 2) {
		INFO("Specify one .prd file");
		return -1;
	}

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	IdTrackHandler<Graph> intIds(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	KmerMapper<K+1, Graph> mapper(g);
	Sequence sequence("");

	long_contigs::LoadFromFile(lc_cfg::get().ds.graph_file, &g, &intIds, sequence, &mapper);

	DataScanner<Graph> dataScanner(g, intIds);
	dataScanner.loadPaired(std::string(argv[1]), pairedIndex);

	PairedInfoChecker checker(g);
	INFO("Symmetic: " << checker.IsSymmetric(pairedIndex));

	return 0;
}


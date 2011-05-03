#include "graphBuilder.hpp"
#include <ext/hash_map>
#include "strobe_reader.hpp"
#include "logging.hpp"
#include "read.hpp"

namespace abruijn {

bool isTrusted(hash_t hash) {
	return true;
}

/**
 * Calculates HTAKE (distinct) minimal hash-values of all k-mers in this read,
 * and puts them into earmarked_hashes
 */
void GraphBuilder::findMinimizers(Sequence s) {
	ha.resize(max(s.size(), ha.size()));
	hashSym.kmers(s, ha);
	for (size_t i = 0; i < htake; i++) {
		hbest[i] = hashing::kMax;
	}
	for (size_t i = 0; i + K <= s.size(); ++i) {
		hash_t hi = ha[i];
		if (!isTrusted(hi)) {
			continue;
		}
		size_t j = htake;
		while (j > 0 && hi < hbest[j - 1]) {
			--j;
		}
		if (j == htake || hi == hbest[j]) {
			continue;
		}
		for (size_t k = htake - 1; k > j; --k) {
			hbest[k] = hbest[k - 1];
		}
		hbest[j] = hi;
	}
	for (size_t i = 0; i < htake && hbest[i] < hashing::kMax; i++) {
		earmarked_hashes.insert(hbest[i]);
	}
}

/**
 * Marks all k-mers in a given read that are locally minimal in
 * a window of size window_size
 */
void GraphBuilder::findLocalMinimizers(Sequence s, size_t window_size) {
	INFO("seq: " << s << " s.size: " << s.size() << " window_size: " << window_size);
	assert(window_size % 2 == 1);

	/// compute hash-values of all the k-mers of a given read
	ha.resize(max(s.size(), ha.size()));
	hashSym.kmers(s, ha);
	INFO("ha.size: " << ha.size());

	/// compute the minimum hash-value in the first window
	assert(window_size <= ha.size());
    hash_t current_min=*(std::min_element(ha.begin(), ha.begin()+window_size));

	for (size_t i = 0; i + window_size < ha.size(); ++i) {
		/// if the current min is (potentially) lost, update it
		if ((i>0 && ha[i-1]==current_min))
			current_min=*(std::min_element(ha.begin()+i, ha.begin()+i+window_size));

		current_min=min(current_min,ha[i+window_size-1]);

		if (ha[i+window_size/2]==current_min)
		{
			earmarked_hashes.insert(current_min);
			INFO(current_min << " is marked");
		}
	}
}


/**
 * If only one k-mer from s is earmarked,
 * this method earmarks 1 more k-mer (with second minimal hash-value)
 */
void GraphBuilder::findSecondMinimizer(Sequence s) {
	hashSym.kmers(s, ha);
	hash_t he = hashing::kMax;
	hash_t hb = hashing::kMax;
	for (size_t i = 0; i + K <= s.size(); ++i) {
		hash_t hi = ha[i];
		if (earmarked_hashes.count(hi)) {
			if (he == hashing::kMax) {
				he = hi;
			} else {
				return;
			}
		} else if (hi < hb) {
			if (!isTrusted(hi)) {
				continue;
			}
			hb = hi;
		}
	}
	assert(hb < hashing::kMax);
	earmarked_hashes.insert(hb);
}

bool Lesser(const Sequence& s) {
	return s < !s;
}

bool LesserK(const Sequence& s, size_t i) {
	return Lesser(s.Subseq(i, i + K));
}

void GraphBuilder::revealTips(Sequence s) {
	hashSym.kmers(s, ha);
	vector<int> index;
	for (size_t i = 0; i + K <= s.size(); ++i) {
		hash_t hi = ha[i];
		if (earmarked_hashes.count(hi)) {
			index.push_back(i);
		}
	}
	for (size_t i = 0; i < index.size(); ++i) {
		hash_t hi = ha[index[i]];
		bool lesser = LesserK(s, index[i]);
		if (i < index.size() - 1) {
			has_right[hi] |= lesser ? 1 : 2;
		}
		if (i > 0) {
			has_right[hi] |= lesser ? 2 : 1;
		}
	}
}

void GraphBuilder::findTipExtensions(Sequence s) {
	hashSym.kmers(s, ha);
	vector<int> index;
	for (size_t i = 0; i + K <= s.size(); ++i) {
		hash_t hi = ha[i];
		if (tips.count(hi)) {
			index.push_back(i);
		}
	}
	for (size_t i = 0; i < index.size(); ++i) {
		bool lesser = LesserK(s, index[i]);
		int l = lesser ? 1 : 2;
		hash_t hi = ha[index[i]];
		TRACE(index[i] << " " << l << " " << (int) tips[hi]);
		size_t low, high;
		if (tips[hi] == l) {
			low = 0;
			high = index[i];
		} else {
			low = index[i] + 1;
			high = s.size() + 1 - K;
		}
		int x = 0;
		for (size_t j = low; j < high; j++) {
			hash_t hj = ha[j];
			if (!isTrusted(hj)) {
				continue;
			}
			assert(!earmarked_hashes.count(hj));
			has_right[hj] = 0;
			tip_extensions[hi].insert(hj);
			x++;
		}
		TRACE(x << " possible continuations");
	}
}

void GraphBuilder::lookRight(Sequence s) {
	hashSym.kmers(s, ha);
	vector<size_t> index;
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (has_right.count(hi)) {
			index.push_back(i);
		}
	}
	if (index.size() == 0) {
		return;
	}
	size_t low = 0;
	size_t high = s.size() - K;
	while (!earmarked_hashes.count(ha[low])) low++;
	while (!earmarked_hashes.count(ha[high])) high--;
	for (size_t i = 0; i < index.size(); ++i) {
		hash_t hi = ha[index[i]];
		bool lesser = LesserK(s, i);
		if (index[i] < high) {
			has_right[hi] |= lesser ? 1 : 2;
		}
		if (index[i] > low) {
			has_right[hi] |= lesser ? 2 : 1;
		}
	}
}

//void extendTip(Sequence s) {
//
//}

void GraphBuilder::addToGraph(Sequence s) {
	vector<abruijn::Vertex*> vs;
	vector<size_t> index;
	hashSym.kmers(s, ha);
	for (size_t i = 0; i + K <= s.size(); i++) {
		if (earmarked_hashes.find(ha[i]) != earmarked_hashes.end()) {
			vs.push_back(graph.getVertex(s.Subseq(i, i + K)));
			index.push_back(i);
		}
	}
	for (size_t i = 0; i + 1 < vs.size(); ++i) {
		graph.addEdge(vs[i], vs[i + 1], s.Subseq(index[i], index[i + 1] + K));
	}
}

void GraphBuilder::build() {
	Read r;

	INFO("===== Finding " << htake << " minimizers in each read... =====");
	srw_.reset();
	for (size_t i = 0; !srw_.eof(); ++i) {
		srw_ >> r;
		if (mode_ & 1) {
			findMinimizers(r.getSequence());
		} else {
			findLocalMinimizers(r.getSequence(), 51);
		}
		VERBOSE(i, " single reads");
	}
	INFO("Done: " << earmarked_hashes.size() << " earmarked hashes");

	if ((mode_ & 2) && (htake == 1)) {
		INFO("===== Finding second minimizers... =====");
		srw_.reset();
		for (size_t i = 0; !srw_.eof(); ++i) {
			srw_ >> r;
			findSecondMinimizer(r.getSequence());
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << earmarked_hashes.size() << " earmarked hashes");
	}

	for(;;) {
		size_t eh = earmarked_hashes.size();
		has_right.clear();
		tips.clear();
		tip_extensions.clear();

		INFO("===== Revealing tips... =====");
		srw_.reset();
		for (size_t i = 0; !srw_.eof(); ++i) {
			srw_ >> r;
			revealTips(r.getSequence());
			VERBOSE(i, " single reads");
		}
		for (map<hash_t, char>::iterator it = has_right.begin(); it != has_right.end(); ++it) {
			if (it->second != 3) {
				TRACE(it->first << " " << (int) it->second);
				tips.insert(*it);
			}
		}
		has_right.clear();
		INFO("Done: " << tips.size() << " tips.");

		INFO("===== Finding tip extensions... =====");
		srw_.reset();
		for (size_t i = 0; !srw_.eof(); ++i) {
			srw_ >> r;
			findTipExtensions(r.getSequence());
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << has_right.size() << " possible tip extensions");

		INFO("===== Looking to the right... =====");
		srw_.reset();
		for (size_t i = 0; !srw_.eof(); ++i) {
			srw_ >> r;
			lookRight(r.getSequence());
			VERBOSE(i, " single reads");
		}
		for (map<hash_t, set<hash_t> >::iterator it = tip_extensions.begin(); it != tip_extensions.end(); ++it) {
			bool ok = false;
			hash_t found = hashing::kMax;
			for (set<hash_t>::iterator ext = it->second.begin(); ext != it->second.end(); ext++) {
				if (earmarked_hashes.count(*ext)) {
					ok = true;
					break;
				}
				if (has_right[*ext] == 3 && found == hashing::kMax) {
					found = *ext;
				}
			}
			if (ok) {
				continue;
			}
			earmarked_hashes.insert(found);
		}
		INFO("Done: " << eh << " -> " << earmarked_hashes.size() << " earmarked hashes");
		if (eh == earmarked_hashes.size()) {
			break;
		}
	}

	INFO("===== Adding reads to graph as paths... =====");
	srw_.reset();
	for (size_t i = 0; !srw_.eof(); ++i) {
		srw_ >> r;
		addToGraph(r.getSequence());
		VERBOSE(i, " single reads");
	}
	INFO("Done: " << graph.vertices.size() << " vertices");

	INFO("===== Condensing-A graph... =====");
	graph.Condense();
	INFO("Done: " << graph.vertices.size() << " vertices");

	return;
}

}

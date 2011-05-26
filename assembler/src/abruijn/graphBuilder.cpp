#include "graphBuilder.hpp"
#include <ext/hash_map>
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
	for (size_t i = 0; i < htake_; i++) {
		hbest[i] = hashing::kMax;
	}
	for (size_t i = 0; i + K <= s.size(); ++i) {
		hash_t hi = ha[i];
		if (!isTrusted(hi)) {
			continue;
		}
		size_t j = htake_;
		while (j > 0 && hi < hbest[j - 1]) {
			--j;
		}
		if (j == htake_ || hi == hbest[j]) {
			continue;
		}
		for (size_t k = htake_ - 1; k > j; --k) {
			hbest[k] = hbest[k - 1];
		}
		hbest[j] = hi;
	}
	for (size_t i = 0; i < htake_ && hbest[i] < hashing::kMax; i++) {
		earmarked_hashes.insert(hbest[i]);
	}
}

/**
 * Marks all k-mers in a given read that are locally minimal in
 * a window of size window_size
 */
void GraphBuilder::findLocalMinimizers(Sequence s, size_t window_size) {
	//INFO("seq: " << s << " s.size: " << s.size() << " window_size: " << window_size);
	assert(window_size % 2 == 1);

	/// compute hash-values of all the k-mers of a given read
	ha.resize(max(s.size(), ha.size()));
	hashSym.kmers(s, ha);
	//INFO("ha.size: " << ha.size());

	/// compute the minimum hash-value in the first window
	window_size=min(window_size, ha.size());
	assert(window_size <= ha.size());
    hash_t current_min=*(std::min_element(ha.begin(), ha.begin()+window_size));

	for (size_t i = 0; i + window_size < ha.size(); ++i) {
		/// if the current min is (potentially) lost, update it
		if ((i>0 && ha[i-1]==current_min))
			current_min=*(std::min_element(ha.begin()+i, ha.begin()+i+window_size));

		current_min=min(current_min,ha[i+window_size-1]);

		if (ha[i+window_size/2]==current_min && isTrusted(current_min))
		{
			earmarked_hashes.insert(current_min);
			//INFO(current_min << " is marked");
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
	// SK: btw, what if a read contains no trusted k-mers?
	// we should ensure that even in this (unlikely) case
	// the code does not break
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
	vector<size_t> index;
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
			size_t dist = (j < index[i]) ? (index[i] - j) : (j - index[i]);
			tip_extensions[hi].insert(make_pair(hj, dist));
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
		bool lesser = LesserK(s, index[i]);
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
			vs.push_back(graph_.getVertex(s.Subseq(i, i + K)));
			index.push_back(i);
		}
	}
	for (size_t i = 0; i + 1 < vs.size(); ++i) {
		graph_.addEdge(vs[i], vs[i + 1], s.Subseq(index[i], index[i + 1] + K));
	}
}

}

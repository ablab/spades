#include "position_kmer.hpp"
#include "kmerno.hpp"

bool KMerNo::equal(const KMerNo & kmerno) const {
	return ( (strncmp( PositionKMer::blob + index, PositionKMer::blob + kmerno.index, K) == 0) );
}

string KMerNo::str() const {
	string res = "";
	for (uint32_t i = 0; i < K; ++i) {
		res += PositionKMer::blob[ index + i ];
	}
	return res;
}

bool KMerNo::less(const KMerNo &r) const {
	return ( strncmp( PositionKMer::blob + index, PositionKMer::blob + r.index, K) < 0 );
}

bool KMerNo::greater(const KMerNo &r) const {
	return ( strncmp( PositionKMer::blob + index, PositionKMer::blob + r.index, K) > 0 );
}

inline char my_dignucl(char c) {
	switch(c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'N': return 4;
		default: assert(false); return -1; // never happens
	}
}

size_t KMerNo::hash::operator() (const KMerNo &kn) const {
	size_t h = 479;
	for (size_t i = 0; i < K; i++) {
		h = ((h << 5) - h) + my_dignucl(PositionKMer::blob[kn.index + i]);
	}
	return h;
}

bool KMerNo::are_equal::operator() (const KMerNo &l, const KMerNo &r) const {
	return ( (strncmp( PositionKMer::blob + l.index, PositionKMer::blob + r.index, K) == 0) );
}


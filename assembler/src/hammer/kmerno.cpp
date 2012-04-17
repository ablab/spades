//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "position_kmer.hpp"
#include "kmerno.hpp"

bool KMerNo::equal(const KMerNo & kmerno) const {
	return ( (strncmp( Globals::blob + index, Globals::blob + kmerno.index, K) == 0) );
}

bool KMerNo::equal(const KMerCount & kmc) const {
	return ( (strncmp( Globals::blob + index, Globals::blob + kmc.first.start(), K) == 0) );
}

string KMerNo::str() const {
	string res = "";
	for (uint32_t i = 0; i < K; ++i) {
		res += Globals::blob[ index + i ];
	}
	return res;
}

bool KMerNo::less(const KMerNo &r) const {
	return ( strncmp( Globals::blob + index, Globals::blob + r.index, K) < 0 );
}

bool KMerNo::greater(const KMerNo &r) const {
	return ( strncmp( Globals::blob + index, Globals::blob + r.index, K) > 0 );
}

inline char my_dignucl(char c) {
	switch(c) {
		case 'A': return 1;
		case 'C': return 2;
		case 'G': return 3;
		case 'T': return 4;
		case 'N': return 5;
		default: assert(false); return -1; // never happens
	}
}

uint64_t KMerNo::new_hash( hint_t index ) {
	hint_t res = 0;
	for (int i = K-1; i >=0; --i) {
		res = ( res*KMERNO_HASH_Q + my_dignucl(Globals::blob[index + i]) ) % KMERNO_HASH_MODULUS;
	}
	return res;
}

uint64_t KMerNo::next_hash( uint64_t old_hash, hint_t new_index ) {
	return (( (old_hash - my_dignucl(Globals::blob[new_index-1]))*KMERNO_HASH_Q_INV
		+ my_dignucl(Globals::blob[new_index+K-1])*KMERNO_HASH_Q_POW_K_MINUS_ONE ) % KMERNO_HASH_MODULUS);
}

uint64_t KMerNo::hash::operator() (const KMerNo &kn) const {
	size_t h = 239;
	for (size_t i = 0; i < K; i++) {
		h = ((h << 5) - h) + Globals::blob[kn.index+i];
	}
	return h;
}

uint64_t KMerNo::string_hash::operator() (const string &kn) const {
	size_t h = 239;
	for (size_t i = 0; i < K; i++) {
		h = ((h << 5) - h) + kn[i];
	}
	return h;
}

bool KMerNo::are_equal::operator() (const KMerNo &l, const KMerNo &r) const {
	return ( (strncmp( Globals::blob + l.index, Globals::blob + r.index, K) == 0) );
}

bool KMerNo::is_less::operator() (const KMerNo &l, const KMerNo &r) const {
	return ( (strncmp( Globals::blob + l.index, Globals::blob + r.index, K) < 0) );
}

bool KMerNo::is_less_kmercount::operator() (const KMerCount &l, const KMerCount &r) const {
	return ( (strncmp( Globals::blob + l.first.start(), Globals::blob + r.first.start(), K) < 0) );
}

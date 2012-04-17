//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "hash.hpp"

namespace hashing {

hash_t power(size_t k) {
	hash_t h = 1;
	for (size_t i = 0; i < k; i++) {
		h = mult(h);
	}
	return h;
}

}

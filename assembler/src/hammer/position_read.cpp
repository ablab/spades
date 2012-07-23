//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "position_kmer.hpp"
#include "position_read.hpp"

char PositionRead::at(size_t pos) const {
	return Globals::blob[start_ + pos];
}

char PositionRead::operator[](size_t pos) const {
	return Globals::blob[start_ + pos];
}

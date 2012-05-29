//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_stat.cpp
 *
 *  Created on: 10.03.2012
 *      Author: snikolenko
 */

#include "standard.hpp"
#include <iostream>
#include <fstream>
#include <queue>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/format.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>

#include <time.h>
#include <sys/resource.h>
#include <iomanip>
#include <algorithm>

#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
#include "subkmers.hpp"
#include "hammer_tools.hpp"

#include "kmer_stat.hpp"


char getQual(const KMerCount & kmc, int i) {
	if (Globals::use_common_quality) return Globals::common_quality * kmc.second.count;
	if (kmc.second.count == 1) {
		// if (Globals::blobquality[kmc.first.start() + i] < Globals::char_offset + 2) cout << "Zero! " << kmc.first.str() << endl;
		return Globals::blobquality[kmc.first.start() + i] - Globals::char_offset;
	}
	else return kmc.second.qual[i];
}


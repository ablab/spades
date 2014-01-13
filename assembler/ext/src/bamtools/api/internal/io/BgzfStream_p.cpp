// ***************************************************************************
// BgzfStream_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 17 January 2012(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#include "bamtools/api/BamAux.h"
#include "bamtools/api/BamConstants.h"
#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include "zlib.h"

#include <cstring>
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;

// ---------------------------
// BgzfStream implementation
// ---------------------------

// checks BGZF block header
bool BgzfStream::CheckBlockHeader(char* header) {
    return (header[0] == Constants::GZIP_ID1 &&
            header[1] == Constants::GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & Constants::FLG_FEXTRA) != 0 &&
            BamTools::UnpackUnsignedShort(&header[10]) == Constants::BGZF_XLEN &&
            header[12] == Constants::BGZF_ID1 &&
            header[13] == Constants::BGZF_ID2 &&
            BamTools::UnpackUnsignedShort(&header[14]) == Constants::BGZF_LEN );
}

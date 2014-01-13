// ***************************************************************************
// ParallelBgzfStream_p.cpp (c) 2011 Derek Barnett
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
#include "api/internal/io/ParallelBgzfStream_p.h"
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
// ParallelBgzfStream implementation
// ---------------------------

// constructor
ParallelBgzfStream::ParallelBgzfStream(void)
  : m_pbgzf(NULL)
  , m_numThreads(0)
{ }

ParallelBgzfStream::ParallelBgzfStream(int32_t numThreads) 
  : m_pbgzf(NULL)
  , m_numThreads(numThreads)
{ }

// destructor
ParallelBgzfStream::~ParallelBgzfStream(void) {
    Close();
}

// closes BGZF file
void ParallelBgzfStream::Close(void) {

    // skip if no device open
    if ( NULL == m_pbgzf) return;

    // close device
    pbgzf_close(m_pbgzf);
    m_pbgzf = NULL;
}

bool ParallelBgzfStream::IsOpen(void) const {
    if ( NULL == m_pbgzf )
        return false;
    return true;
}

void ParallelBgzfStream::Open(const string& filename, const IBamIODevice::OpenMode mode) {
    char mode_cstr[16] = "\0";
    // close current device if necessary
    Close();
    BT_ASSERT_X( (NULL == m_pbgzf), "PgzfStream::Open() - unable to properly close previous IO device" );

    switch (mode) {
        case IBamIODevice::ReadOnly:
            strcat(mode_cstr, "r");
            break;
        case IBamIODevice::WriteOnly:
            strcat(mode_cstr, "w");
            if(m_isWriteCompressed) 
                strcat(mode_cstr, "b");
            break;
        default:
            // TODO: why would this every happen?  Don't call me this way.
            throw BamException("PgzfStream::Open", "could not open BGZF stream");
        break;
    }
            
    pbgzf_set_num_threads_per(m_numThreads); // TODO: this needs to be set on a per "m_pbgzf" basis
    m_pbgzf = pbgzf_open(filename.c_str(), mode_cstr); 
    // TODO: asserts

}

// reads BGZF data into a byte buffer
size_t ParallelBgzfStream::Read(char* data, const size_t dataLength) {

    if ( dataLength == 0 )
        return 0;

    // TODO: asserts
    return pbgzf_read(m_pbgzf, (void*)data, (int)dataLength);
}

// seek to position in BGZF file
void ParallelBgzfStream::Seek(const int64_t& position) {

    // TODO: asserts
    pbgzf_seek(m_pbgzf, position, SEEK_SET);
}

// get file position in BGZF file
int64_t ParallelBgzfStream::Tell(void) const {
    if ( NULL == m_pbgzf )
        return 0;
    return pbgzf_tell(m_pbgzf);
}

// writes the supplied data into the BGZF buffer
size_t ParallelBgzfStream::Write(const char* data, const size_t dataLength) {
    
    // TODO: asserts
    return (size_t)pbgzf_write(m_pbgzf, data, dataLength);
}

// ***************************************************************************
// ParallelBgzfStream_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 17 January 2012(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#ifndef PARALLELBGZFSTREAM_P_H
#define PARALLELBGZFSTREAM_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "bamtools/api/api_global.h"
#include "bamtools/api/BamAux.h"
#include "bamtools/api/IBamIODevice.h"
#include "api/internal/io/pbgzf/pbgzf.h"
#include "api/internal/io/BgzfStream_p.h"
#include <string>

namespace BamTools {
namespace Internal {

class ParallelBgzfStream : public BgzfStream {

    // constructor & destructor
    public:
        ParallelBgzfStream(void);
        ParallelBgzfStream(int numThreads);
        ~ParallelBgzfStream(void);

    // main interface methods
    public:
        // closes BGZF file
        void Close(void);
        // returns true if ParallelBgzfStream open for IO
        bool IsOpen(void) const;
        // opens the BGZF file
        void Open(const std::string& filename, const IBamIODevice::OpenMode mode);
        // reads BGZF data into a byte buffer
        size_t Read(char* data, const size_t dataLength);
        // seek to position in BGZF file
        void Seek(const int64_t& position);
        // enable/disable compressed output
        void SetWriteCompressed(bool ok);
        // get file position in BGZF file
        int64_t Tell(void) const;
        // writes the supplied data into the BGZF buffer
        size_t Write(const char* data, const size_t dataLength);

    // data members
    public:

        PBGZF *m_pbgzf;
        int32_t m_numThreads;
};

} // namespace Internal
} // namespace BamTools

#endif // PARALLELBGZFSTREAM_P_H

// ***************************************************************************
// BgzfStream_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 17 January 2012(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#ifndef BGZFSTREAM_P_H
#define BGZFSTREAM_P_H

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
#include <string>

namespace BamTools {
namespace Internal {

// abstract class
class BgzfStream {

    // main interface methods
    public:
        // closes BGZF file
        virtual void Close(void) = 0;
        // returns true if BgzfStream open for IO
        virtual bool IsOpen(void) const = 0;
        // opens the BGZF file
        virtual void Open(const std::string& filename, const IBamIODevice::OpenMode mode) = 0;
        // reads BGZF data into a byte buffer
        virtual size_t Read(char* data, const size_t dataLength) = 0;
        // seek to position in BGZF file
        virtual void Seek(const int64_t& position) = 0;
        // sets IO device (closes previous, if any, but does not attempt to open)
        //void SetIODevice(IBamIODevice* device) = { return; };
        // enable/disable compressed output
        void SetWriteCompressed(bool ok) { m_isWriteCompressed = ok; };
        // get file position in BGZF file
        virtual int64_t Tell(void) const = 0;
        // writes the supplied data into the BGZF buffer
        virtual size_t Write(const char* data, const size_t dataLength) = 0;


    // static 'utility' methods
    public:
        // checks BGZF block header
        static bool CheckBlockHeader(char* header);

    protected:
        bool m_isWriteCompressed;
};

} // namespace Internal
} // namespace BamTools

#endif // BGZFSTREAM_P_H

// ***************************************************************************
// BamReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 18 November 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_P_H
#define BAMREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "bamtools/api/BamAlignment.h"
#include "bamtools/api/BamIndex.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamHeader.h"
#include "api/internal/bam/BamHeader_p.h"
#include "api/internal/bam/BamRandomAccessController_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include <string>

namespace BamTools {
namespace Internal {

class BamReaderPrivate {

    // ctor & dtor
    public:
        BamReaderPrivate(BamReader* parent);
        ~BamReaderPrivate(void);

    // BamReader interface
    public:

        // file operations
        bool Close(void);
        const std::string Filename(void) const;
        bool IsOpen(void) const;
        bool Open(const std::string& filename);
        bool Rewind(void);
        bool SetRegion(const BamRegion& region);

        // access alignment data
        bool GetNextAlignment(BamAlignment& alignment);
        bool GetNextAlignmentCore(BamAlignment& alignment);

        // access auxiliary data
        std::string GetHeaderText(void) const;
        const SamHeader& GetConstSamHeader(void) const;
        SamHeader GetSamHeader(void) const;
        int GetReferenceCount(void) const;
        const RefVector& GetReferenceData(void) const;
        int GetReferenceID(const std::string& refName) const;

        // index operations
        bool CreateIndex(const BamIndex::IndexType& type);
        bool HasIndex(void) const;
        bool LocateIndex(const BamIndex::IndexType& preferredType);
        bool OpenIndex(const std::string& indexFilename);
        void SetIndex(BamIndex* index);

        // error handling
        std::string GetErrorString(void) const;
        void SetErrorString(const std::string& where, const std::string& what);

        void SetParallel(int32_t);

    // internal methods, but available as a BamReaderPrivate 'interface'
    //
    // these methods should only be used by BamTools::Internal classes
    // (currently only used by the BamIndex subclasses)
    public:
        // retrieves header text from BAM file
        void LoadHeaderData(void);
        // retrieves BAM alignment under file pointer
        // (does no overlap checking or character data parsing)
        bool LoadNextAlignment(BamAlignment& alignment);
        // builds reference data structure from BAM file
        bool LoadReferenceData(void);
        // seek reader to file position
        bool Seek(const int64_t& position);
        // return reader's file position
        int64_t Tell(void) const;

    // data members
    public:

        // general BAM file data
        int64_t     m_alignmentsBeginOffset;
        std::string m_filename;
        RefVector   m_references;

        // system data
        bool m_isBigEndian;

        // parent BamReader
        BamReader* m_parent;

        // BamReaderPrivate components
        BamHeader m_header;
        BamRandomAccessController m_randomAccessController;
        BgzfStream *m_stream;

        // error handling
        std::string m_errorString;

        int32_t m_numThreads;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMREADER_P_H

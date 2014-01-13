// ***************************************************************************
// SerialBgzfStream_p.cpp (c) 2011 Derek Barnett
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
#include "api/internal/io/SerialBgzfStream_p.h"
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
// SerialBgzfStream implementation
// ---------------------------

// constructor
SerialBgzfStream::SerialBgzfStream(void)
  : m_blockLength(0)
  , m_blockOffset(0)
  , m_blockAddress(0)
  , m_isWriteCompressed(true)
  , m_device(0)
  , m_uncompressedBlock(Constants::BGZF_DEFAULT_BLOCK_SIZE)
  , m_compressedBlock(Constants::BGZF_MAX_BLOCK_SIZE)
{ }

// destructor
SerialBgzfStream::~SerialBgzfStream(void) {
    Close();
}

// closes BGZF file
void SerialBgzfStream::Close(void) {

    // skip if no device open
    if ( m_device == 0 ) return;

    // if writing to file, flush the current BGZF block,
    // then write an empty block (as EOF marker)
    if ( m_device->IsOpen() && (m_device->Mode() == IBamIODevice::WriteOnly) ) {
        FlushBlock();
        const size_t blockLength = DeflateBlock(0);
        m_device->Write(m_compressedBlock.Buffer, blockLength);
    }

    // close device
    m_device->Close();
    delete m_device;
    m_device = 0;

    // ensure our buffers are cleared out
    m_uncompressedBlock.Clear();
    m_compressedBlock.Clear();

    // reset state
    m_blockLength = 0;
    m_blockOffset = 0;
    m_blockAddress = 0;
    m_isWriteCompressed = true;
}

// compresses the current block
size_t SerialBgzfStream::DeflateBlock(int32_t blockLength) {

    // initialize the gzip header
    char* buffer = m_compressedBlock.Buffer;
    memset(buffer, 0, 18);
    buffer[0]  = Constants::GZIP_ID1;
    buffer[1]  = Constants::GZIP_ID2;
    buffer[2]  = Constants::CM_DEFLATE;
    buffer[3]  = Constants::FLG_FEXTRA;
    buffer[9]  = Constants::OS_UNKNOWN;
    buffer[10] = Constants::BGZF_XLEN;
    buffer[12] = Constants::BGZF_ID1;
    buffer[13] = Constants::BGZF_ID2;
    buffer[14] = Constants::BGZF_LEN;

    // set compression level
    const int compressionLevel = ( m_isWriteCompressed ? Z_DEFAULT_COMPRESSION : 0 );

    // loop to retry for blocks that do not compress enough
    int inputLength = blockLength;
    size_t compressedLength = 0;
    const unsigned int bufferSize = Constants::BGZF_MAX_BLOCK_SIZE;

    while ( true ) {

        // initialize zstream values
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)m_uncompressedBlock.Buffer;
        zs.avail_in  = inputLength;
        zs.next_out  = (Bytef*)&buffer[Constants::BGZF_BLOCK_HEADER_LENGTH];
        zs.avail_out = bufferSize -
                       Constants::BGZF_BLOCK_HEADER_LENGTH -
                       Constants::BGZF_BLOCK_FOOTER_LENGTH;

        // initialize the zlib compression algorithm
        int status = deflateInit2(&zs,
                                  compressionLevel,
                                  Z_DEFLATED,
                                  Constants::GZIP_WINDOW_BITS,
                                  Constants::Z_DEFAULT_MEM_LEVEL,
                                  Z_DEFAULT_STRATEGY);
        if ( status != Z_OK )
            throw BamException("SerialBgzfStream::DeflateBlock", "zlib deflateInit2 failed");

        // compress the data
        status = deflate(&zs, Z_FINISH);

        // if not at stream end
        if ( status != Z_STREAM_END ) {

            deflateEnd(&zs);

            // there was not enough space available in buffer
            // try to reduce the input length & re-start loop
            if ( status == Z_OK ) {
                inputLength -= 1024;
                if ( inputLength < 0 )
                    throw BamException("SerialBgzfStream::DeflateBlock", "input reduction failed");
                continue;
            }

            throw BamException("SerialBgzfStream::DeflateBlock", "zlib deflate failed");
        }

        // finalize the compression routine
        status = deflateEnd(&zs);
        if ( status != Z_OK )
            throw BamException("SerialBgzfStream::DeflateBlock", "zlib deflateEnd failed");

        // update compressedLength
        compressedLength = zs.total_out +
                           Constants::BGZF_BLOCK_HEADER_LENGTH +
                           Constants::BGZF_BLOCK_FOOTER_LENGTH;
        if ( compressedLength > Constants::BGZF_MAX_BLOCK_SIZE )
            throw BamException("SerialBgzfStream::DeflateBlock", "deflate overflow");

        // quit while loop
        break;
    }

    // store the compressed length
    BamTools::PackUnsignedShort(&buffer[16], static_cast<uint16_t>(compressedLength - 1));

    // store the CRC32 checksum
    uint32_t crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)m_uncompressedBlock.Buffer, inputLength);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 8], crc);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

    // ensure that we have less than a block of data left
    int remaining = blockLength - inputLength;
    if ( remaining > 0 ) {
        if ( remaining > inputLength )
            throw BamException("SerialBgzfStream::DeflateBlock", "after deflate, remainder too large");
        memcpy(m_uncompressedBlock.Buffer, m_uncompressedBlock.Buffer + inputLength, remaining);
    }

    // update block data
    m_blockOffset = remaining;

    // return result
    return compressedLength;
}

// flushes the data in the BGZF block
void SerialBgzfStream::FlushBlock(void) {

    BT_ASSERT_X( m_device, "SerialBgzfStream::FlushBlock() - attempting to flush to null device" );

    // flush all of the remaining blocks
    while ( m_blockOffset > 0 ) {

        // compress the data block
        const size_t blockLength = DeflateBlock(m_blockOffset);

        // flush the data to our output device
        const int64_t numBytesWritten = m_device->Write(m_compressedBlock.Buffer, blockLength);

        // check for device error
        if ( numBytesWritten < 0 ) {
            const string message = string("device error: ") + m_device->GetErrorString();
            throw BamException("SerialBgzfStream::FlushBlock", message);
        }

        // check that we wrote expected numBytes
        if ( numBytesWritten != static_cast<int64_t>(blockLength) ) {
            stringstream s("");
            s << "expected to write " << blockLength
              << " bytes during flushing, but wrote " << numBytesWritten;
            throw BamException("SerialBgzfStream::FlushBlock", s.str());
        }

        // update block data
        m_blockAddress += blockLength;
    }
}

// decompresses the current block
size_t SerialBgzfStream::InflateBlock(const size_t& blockLength) {

    // setup zlib stream object
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)m_compressedBlock.Buffer + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)m_uncompressedBlock.Buffer;
    zs.avail_out = Constants::BGZF_DEFAULT_BLOCK_SIZE;

    // initialize
    int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);
    if ( status != Z_OK )
        throw BamException("SerialBgzfStream::InflateBlock", "zlib inflateInit failed");

    // decompress
    status = inflate(&zs, Z_FINISH);
    if ( status != Z_STREAM_END ) {
        inflateEnd(&zs);
        throw BamException("SerialBgzfStream::InflateBlock", "zlib inflate failed");
    }

    // finalize
    status = inflateEnd(&zs);
    if ( status != Z_OK ) {
        inflateEnd(&zs);
        throw BamException("SerialBgzfStream::InflateBlock", "zlib inflateEnd failed");
    }

    // return result
    return zs.total_out;
}

bool SerialBgzfStream::IsOpen(void) const {
    if ( m_device == 0 )
        return false;
    return m_device->IsOpen();
}

void SerialBgzfStream::Open(const string& filename, const IBamIODevice::OpenMode mode) {

    // close current device if necessary
    Close();
    BT_ASSERT_X( (m_device == 0), "SerialBgzfStream::Open() - unable to properly close previous IO device" );

    // retrieve new IO device depending on filename
    m_device = BamDeviceFactory::CreateDevice(filename);
    BT_ASSERT_X( m_device, "SerialBgzfStream::Open() - unable to create IO device from filename" );

    // if device fails to open
    if ( !m_device->Open(mode) ) {
        const string deviceError = m_device->GetErrorString();
        const string message = string("could not open BGZF stream: \n\t") + deviceError;
        throw BamException("SerialBgzfStream::Open", message);
    }
}

// reads BGZF data into a byte buffer
size_t SerialBgzfStream::Read(char* data, const size_t dataLength) {

    if ( dataLength == 0 )
        return 0;

    // if stream not open for reading
    BT_ASSERT_X( m_device, "SerialBgzfStream::Read() - trying to read from null device");
    if ( !m_device->IsOpen() || (m_device->Mode() != IBamIODevice::ReadOnly) )
        return 0;

    // read blocks as needed until desired data length is retrieved
    char* output = data;
    size_t numBytesRead = 0;
    while ( numBytesRead < dataLength ) {

        // determine bytes available in current block
        int bytesAvailable = m_blockLength - m_blockOffset;

        // read (and decompress) next block if needed
        if ( bytesAvailable <= 0 ) {
            ReadBlock();
            bytesAvailable = m_blockLength - m_blockOffset;
            if ( bytesAvailable <= 0 )
                break;
        }

        // copy data from uncompressed source buffer into data destination buffer
        const size_t copyLength = min( (dataLength-numBytesRead), (size_t)bytesAvailable );
        memcpy(output, m_uncompressedBlock.Buffer + m_blockOffset, copyLength);

        // update counters
        m_blockOffset += copyLength;
        output        += copyLength;
        numBytesRead  += copyLength;
    }

    // update block data
    if ( m_blockOffset == m_blockLength ) {
        m_blockAddress = m_device->Tell();
        m_blockOffset  = 0;
        m_blockLength  = 0;
    }

    // return actual number of bytes read
    return numBytesRead;
}

// reads a BGZF block
void SerialBgzfStream::ReadBlock(void) {

    BT_ASSERT_X( m_device, "SerialBgzfStream::ReadBlock() - trying to read from null IO device");

    // store block's starting address
    const int64_t blockAddress = m_device->Tell();

    // read block header from file
    char header[Constants::BGZF_BLOCK_HEADER_LENGTH];
    int64_t numBytesRead = m_device->Read(header, Constants::BGZF_BLOCK_HEADER_LENGTH);

    // check for device error
    if ( numBytesRead < 0 ) {
        const string message = string("device error: ") + m_device->GetErrorString();
        throw BamException("SerialBgzfStream::ReadBlock", message);
    }

    // if block header empty
    if ( numBytesRead == 0 ) {
        m_blockLength = 0;
        return;
    }

    // if block header invalid size
    if ( numBytesRead != static_cast<int8_t>(Constants::BGZF_BLOCK_HEADER_LENGTH) )
        throw BamException("SerialBgzfStream::ReadBlock", "invalid block header size");

    // validate block header contents
    if ( !BgzfStream::CheckBlockHeader(header) )
        throw BamException("SerialBgzfStream::ReadBlock", "invalid block header contents");

    // copy header contents to compressed buffer
    const size_t blockLength = BamTools::UnpackUnsignedShort(&header[16]) + 1;
    memcpy(m_compressedBlock.Buffer, header, Constants::BGZF_BLOCK_HEADER_LENGTH);

    // read remainder of block
    const size_t remaining = blockLength - Constants::BGZF_BLOCK_HEADER_LENGTH;
    numBytesRead = m_device->Read(&m_compressedBlock.Buffer[Constants::BGZF_BLOCK_HEADER_LENGTH], remaining);

    // check for device error
    if ( numBytesRead < 0 ) {
        const string message = string("device error: ") + m_device->GetErrorString();
        throw BamException("SerialBgzfStream::ReadBlock", message);
    }

    // check that we read in expected numBytes
    if ( numBytesRead != static_cast<int64_t>(remaining) )
        throw BamException("SerialBgzfStream::ReadBlock", "could not read data from block");

    // decompress block data
    const size_t newBlockLength = InflateBlock(blockLength);

    // update block data
    if ( m_blockLength != 0 )
        m_blockOffset = 0;
    m_blockAddress = blockAddress;
    m_blockLength  = newBlockLength;
}

// seek to position in BGZF file
void SerialBgzfStream::Seek(const int64_t& position) {

    BT_ASSERT_X( m_device, "SerialBgzfStream::Seek() - trying to seek on null IO device");

    // skip if device is not open
    if ( !IsOpen() ) return;

    // determine adjusted offset & address
    int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;

    // attempt seek in file
    if ( m_device->IsRandomAccess() && m_device->Seek(blockAddress) ) {

        // update block data & return success
        m_blockLength  = 0;
        m_blockAddress = blockAddress;
        m_blockOffset  = blockOffset;
    }
    else {
        stringstream s("");
        s << "unable to seek to position: " << position;
        throw BamException("SerialBgzfStream::Seek", s.str());
    }
}

// get file position in BGZF file
int64_t SerialBgzfStream::Tell(void) const {
    if ( !IsOpen() )
        return 0;
    return ( (m_blockAddress << 16) | (m_blockOffset & 0xFFFF) );
}

// writes the supplied data into the BGZF buffer
size_t SerialBgzfStream::Write(const char* data, const size_t dataLength) {

    BT_ASSERT_X( m_device, "SerialBgzfStream::Write() - trying to write to null IO device");
    BT_ASSERT_X( (m_device->Mode() == IBamIODevice::WriteOnly),
                 "SerialBgzfStream::Write() - trying to write to non-writable IO device");

    // skip if file not open for writing
    if ( !IsOpen() )
        return 0;

    // write blocks as needed til all data is written
    size_t numBytesWritten = 0;
    const char* input = data;
    const size_t blockLength = Constants::BGZF_DEFAULT_BLOCK_SIZE;
    while ( numBytesWritten < dataLength ) {

        // copy data contents to uncompressed output buffer
        unsigned int copyLength = min(blockLength - m_blockOffset, dataLength - numBytesWritten);
        char* buffer = m_uncompressedBlock.Buffer;
        memcpy(buffer + m_blockOffset, input, copyLength);

        // update counter
        m_blockOffset   += copyLength;
        input           += copyLength;
        numBytesWritten += copyLength;

        // flush (& compress) output buffer when full
        if ( m_blockOffset == static_cast<int32_t>(blockLength) )
            FlushBlock();
    }

    // return actual number of bytes written
    return numBytesWritten;
}

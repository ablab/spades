/*
Copyright 2011 Convey Computer Corporation (info@conveycomputer.com)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "globals.h"
#include "tightString.h"
#include "readSet.h"
#include "binarySequences.h"
#include "utility.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include "../third-party/zlib-1.2.3/Win32/include/zlib.h"
#else
#include "../third-party/zlib-1.2.3/zlib.h"
#endif

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

// write defines, typedefs, and protos
#define WRITE_BUF_SHFT		16	// byte shift and mask
#define WRITE_BUF_SIZE		(1<<WRITE_BUF_SHFT)
#define WRITE_BUF_MASK		(WRITE_BUF_SIZE-1)
#define SHORT_NUCL_LENGTH	128	// Nucleotide length (2-bits each)

void computeSecondInPair(ReadSet * reads);

FILE *openCnySeqForRead(const char *fileName, CnyUnifiedSeqFileHeader *seqFileHeader)
{
	FILE *pFile;
	if ((pFile = fopen(fileName, "rb")) == 0) {
		velvetLog("Unable to open %s for reading\n", fileName);
		return NULL;
	}

	if (fread(seqFileHeader, sizeof(*seqFileHeader), 1, pFile) != 1) {
		velvetLog("Unable to read file %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (strncmp((char *)&seqFileHeader->m_magic, "CSQ0", 4) != 0) {
		velvetLog("Unknown format for file %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (seqFileHeader->m_bFileWriteCompleted == false) {
		velvetLog("Corrupted file, %s\n", fileName);
		fclose(pFile);
		return NULL;
	}

	if (seqFileHeader->m_numCategories > CATEGORIES) {
		velvetLog("File %s has %d categories, please rebuild velvet to match\n", fileName, seqFileHeader->m_numCategories);
		fclose(pFile);
		return NULL;
	}

#ifdef COLOR
	if (!seqFileHeader->m_bColor) {
		velvetLog("File %s does not specify color, please rebuild velvet to match\n", fileName);
		fclose(pFile);
		return NULL;
	}
#else
	if (seqFileHeader->m_bColor) {
		velvetLog("File %s specifies color, please rebuild velvet to match\n", fileName);
		fclose(pFile);
		return NULL;
	}
#endif
	return pFile;
}

static boolean refillCnySeqReadBuffer(SequencesReader *seqReadInfo) 
{
	uint64_t readLen = (USF_READ_BUF_SIZE < seqReadInfo->m_unifiedSeqFileHeader.m_seqNuclStoreSize - seqReadInfo->m_readBufPos) ?
		USF_READ_BUF_SIZE : seqReadInfo->m_unifiedSeqFileHeader.m_seqNuclStoreSize - seqReadInfo->m_readBufPos;

	if (readLen == 0)
		return false;

	if (fread(seqReadInfo->m_pReadBuffer, (uint32_t)readLen, 1, seqReadInfo->m_pFile) != 1) {
		velvetLog("Unable to read file\n");
		exit (1);
	}

	seqReadInfo->m_pCurrentReadPtr = seqReadInfo->m_pReadBuffer;
	seqReadInfo->m_pReadBufEnd = seqReadInfo->m_pReadBuffer + readLen;
	seqReadInfo->m_readBufPos += readLen;

	if (seqReadInfo->m_pNextReadPtr >= seqReadInfo->m_pReadBufEnd) {
		seqReadInfo->m_pNextReadPtr -= USF_READ_BUF_SIZE;
	}

	return true;
}

static int32_t readCnySeqUint8(SequencesReader *seqReadInfo)
{
	if (seqReadInfo->m_pCurrentReadPtr == seqReadInfo->m_pReadBufEnd && !refillCnySeqReadBuffer(seqReadInfo))
	{
		return -1;
	}

	// printf("m_pCurrentReadPtr %llx\n", (long long) pReadInfo->m_pCurrentReadPtr);
	return *seqReadInfo->m_pCurrentReadPtr++;
}

uint32_t readCnySeqUint32(SequencesReader *seqReadInfo) 
{
	uint32_t data;
	data = 0;
	int i;
	for (i = 0; i < 4; i += 1)
		data |= readCnySeqUint8(seqReadInfo) << (i*8);
	return data;
}

boolean advanceCnySeqCurrentRead(SequencesReader *seqReadInfo)
{
	// Perform consistency check, unused bits of previous sequence should have a fixed pattern
	uint32_t finalNuclOffset = 1;
	if (seqReadInfo->m_bIsRef) {
		finalNuclOffset += (sizeof(seqReadInfo->m_refCnt));
		finalNuclOffset += (sizeof(RefInfo) * seqReadInfo->m_refCnt);
	}
	if ((seqReadInfo->m_currentReadLength & 3) != 0 && ((seqReadInfo->m_pNextReadPtr - finalNuclOffset) >= seqReadInfo->m_pReadBuffer)) {
		uint8_t mask = 0xFF << (seqReadInfo->m_currentReadLength & 3) * 2;
		if ((*(seqReadInfo->m_pNextReadPtr - finalNuclOffset) & mask) != (0xAA & mask)) {
			velvetLog("Cny seq consistency check failed in advance\n");
#ifdef DEBUG
			abort();
#endif
			exit(1);
		}
	}

	seqReadInfo->m_pCurrentReadPtr = seqReadInfo->m_pNextReadPtr;
	seqReadInfo->m_currentNuclReadIdx = 0;

	// clear ref flag before each code check
	seqReadInfo->m_bIsRef = false;
	seqReadInfo->m_refCnt = 0;

	for(;;) {
		int32_t	code = readCnySeqUint8(seqReadInfo);
		// printf("checking code %d\n", code);
		switch (code & 0xc0) {
			case 0x00:	// short sequence
			case 0x40:
				seqReadInfo->m_currentReadLength = code & 0x7f;
				seqReadInfo->m_pNextReadPtr = seqReadInfo->m_pCurrentReadPtr + ((seqReadInfo->m_currentReadLength + 3) >> 2);
				break;
			case 0x80:	// long sequence
				seqReadInfo->m_currentReadLength = readCnySeqUint32(seqReadInfo);
				seqReadInfo->m_pNextReadPtr = seqReadInfo->m_pCurrentReadPtr + ((seqReadInfo->m_currentReadLength + 3) >> 2);
				if (code & 0x20) {
				    // ref info present
				    seqReadInfo->m_bIsRef = true;
				    seqReadInfo->m_pNextReadPtr += (sizeof(seqReadInfo->m_refCnt));
				    // length is updated once count is read
				}
				break;
			case 0xc0:	// new file / category
				if (code == EOF) {
					return false;
				}
				seqReadInfo->m_currCategory = (Category) readCnySeqUint32(seqReadInfo);
				if (seqReadInfo->m_currCategory < 0 || seqReadInfo->m_currCategory > REFERENCE) {
					velvetLog("Illegal category %d\n", (int32_t) seqReadInfo->m_currCategory);
					exit(1);
				}
				continue;
		}

		if (seqReadInfo->m_currentReadLength > seqReadInfo->m_maxSeqLen ||
				seqReadInfo->m_currentReadLength < seqReadInfo->m_minSeqLen) {
			velvetLog("Cny seq consistency check failed, len mismatch\n");
#ifdef DEBUG
			abort();
#endif
			exit(1);
		}

		return true;

	}
}

void resetCnySeqCurrentRead(SequencesReader *seqReadInfo) 
{
	seqReadInfo->m_pReadBufEnd = seqReadInfo->m_pReadBuffer;
	seqReadInfo->m_pNextReadPtr = seqReadInfo->m_pReadBuffer;
	seqReadInfo->m_pCurrentReadPtr = seqReadInfo->m_pReadBuffer;
	seqReadInfo->m_currentReadLength = 0;
	seqReadInfo->m_readBufPos = 0;

	if (fseek(seqReadInfo->m_pFile, sizeof(CnyUnifiedSeqFileHeader), SEEK_SET) < 0) {
		perror("Unable to seek\n");
		exit(1);
	}

	advanceCnySeqCurrentRead(seqReadInfo);
}

void getCnySeqNucl(SequencesReader *seqReadInfo, uint8_t *sequence) {
	uint32_t nuclIdx;
	for (nuclIdx = 0; nuclIdx < seqReadInfo->m_currentReadLength; nuclIdx += 4) {
		sequence[nuclIdx / 4] = (uint8_t)readCnySeqUint8(seqReadInfo);
	}
}

ReadSet *importCnyReadSet(char *filename)
{
	IDnum sequenceCount, sequenceIndex;
	ReadSet *reads;
	uint8_t *tmp;
	Coordinate totalLength = 0;
	int arrayLength;
	SequencesReader seqReadInfo;
	memset(&seqReadInfo, 0, sizeof(seqReadInfo));

	seqReadInfo.m_pFile = openCnySeqForRead(filename, &seqReadInfo.m_unifiedSeqFileHeader);
	seqReadInfo.m_numCategories = seqReadInfo.m_unifiedSeqFileHeader.m_numCategories;
	seqReadInfo.m_minSeqLen = seqReadInfo.m_unifiedSeqFileHeader.m_minSeqLen;
	seqReadInfo.m_maxSeqLen = seqReadInfo.m_unifiedSeqFileHeader.m_maxSeqLen;
	seqReadInfo.m_bIsRef = false;

	if (seqReadInfo.m_pFile != NULL)
		velvetLog("Reading CNY read set file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

//	readInfo.m_pReadBuffer = mallocOrExit(USF_READ_BUF_SIZE, sizeof(*readInfo.m_pReadBuffer));
	seqReadInfo.m_pReadBuffer = mallocOrExit(USF_READ_BUF_SIZE, uint8_t );
	seqReadInfo.m_pCurrentReadPtr = seqReadInfo.m_pReadBufEnd = 0;

	reads = newReadSet();

	resetCnySeqCurrentRead(&seqReadInfo);
	sequenceCount = seqReadInfo.m_unifiedSeqFileHeader.m_sequenceCnt;

	velvetLog("%li sequences found\n", (long) sequenceCount);

	reads->readCount = sequenceCount;

	if (reads->readCount == 0) {
		reads->sequences = NULL;
		reads->categories = NULL;
		free(seqReadInfo.m_pReadBuffer);
		return reads;	
	}

	reads->sequences = NULL;
	reads->categories = callocOrExit(sequenceCount, Category);
	reads->tSequences = mallocOrExit(sequenceCount, TightString);
	// note there is some overhead with the seq store
	reads->tSeqMem = callocOrExit (seqReadInfo.m_unifiedSeqFileHeader.m_seqNuclStoreSize, char);
	tmp = (uint8_t *) reads->tSeqMem;
	uint8_t * arrayEnd = tmp + seqReadInfo.m_unifiedSeqFileHeader.m_seqNuclStoreSize;
	// read all sequence and category info in one pass
	for (sequenceIndex = 0; sequenceIndex < sequenceCount; sequenceIndex += 1) {
		reads->categories[sequenceIndex] = seqReadInfo.m_currCategory;
		if (sizeof(ShortLength) == sizeof(int16_t) && seqReadInfo.m_currentReadLength > SHRT_MAX) {
			velvetLog("Read %li of length %lli, longer than limit %i\n",
					(long) sequenceIndex + 1, (long long) seqReadInfo.m_currentReadLength, SHRT_MAX);
			velvetLog("You should recompile Velvet with the LONGSEQUENCES option.\n");
			exit(1);
		}
		// only use tString to reduce memory use
		reads->tSequences[sequenceIndex].length = seqReadInfo.m_currentReadLength;
		arrayLength = (reads->tSequences[sequenceIndex].length + 3) / 4;
		if ((tmp + arrayLength) > arrayEnd) {
			velvetLog("array location 0x%lx for seq %ld beyond end 0x%lx\n", (uint64_t) tmp, (uint64_t) sequenceIndex, (uint64_t) arrayEnd);
			exit(1);
		}
		totalLength += arrayLength;
		reads->tSequences[sequenceIndex].sequence = tmp;
		getCnySeqNucl(&seqReadInfo, tmp);
		if (seqReadInfo.m_bIsRef) {
			seqReadInfo.m_refCnt = readCnySeqUint32(&seqReadInfo);
			// now the next ptr is advanced
			seqReadInfo.m_pNextReadPtr += (sizeof(RefInfo) * seqReadInfo.m_refCnt);
			RefInfo refElem;
			uint32_t refIdx;
			for (refIdx = 0; refIdx < seqReadInfo.m_refCnt; refIdx++) {
				// not actually used so just read past refs
				refElem.m_referenceID = readCnySeqUint32(&seqReadInfo);
				refElem.m_pos = readCnySeqUint32(&seqReadInfo);
			}
		}
		tmp += arrayLength;
		if (sequenceIndex < sequenceCount) {
			advanceCnySeqCurrentRead(&seqReadInfo);
		}
	}

	fclose(seqReadInfo.m_pFile);
	computeSecondInPair(reads);

	free(seqReadInfo.m_pReadBuffer);
	velvetLog("Done\n");
	return reads;

}

// write routines
#define ADENINE		0
#define CYTOSINE	1
#define GUANINE		2
#define THYMINE		3
#define INVALID         5

static void cnySeqHostBufferFull(SequencesWriter *seqWriteInfo)
{
    // The current Host buffer is full
    switch (seqWriteInfo->m_hostBuffersInUse) {
	case 1:	// buf[0] is full, start using buf[1]
	    seqWriteInfo->m_pHostBufPtr = seqWriteInfo->m_pWriteBuffer[1];
	    seqWriteInfo->m_pHostBufPtrMax = seqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    seqWriteInfo->m_hostBufferFilePos[1] = seqWriteInfo->m_hostBufferFilePos[0] + WRITE_BUF_SIZE;
	    seqWriteInfo->m_hostBuffersInUse = 2;
	    break;
	case 2: // buf[0] and buf[1] are full, start using buf[2]
	    seqWriteInfo->m_pHostBufPtr = seqWriteInfo->m_pWriteBuffer[2];
	    seqWriteInfo->m_pHostBufPtrMax = seqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    seqWriteInfo->m_hostBufferFilePos[2] = seqWriteInfo->m_hostBufferFilePos[1] + WRITE_BUF_SIZE;
	    seqWriteInfo->m_hostBuffersInUse = 3;
	    break;
	case 3: // all three buffers are full, write out buf[2] and reuse
	    if (fseek(seqWriteInfo->m_pFile, seqWriteInfo->m_hostBufferFilePos[2], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }

	    if (fwrite(seqWriteInfo->m_pWriteBuffer[2], WRITE_BUF_SIZE, 1, seqWriteInfo->m_pFile) != 1) {
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    seqWriteInfo->m_pHostBufPtr = seqWriteInfo->m_pWriteBuffer[2];
	    seqWriteInfo->m_pHostBufPtrMax = seqWriteInfo->m_pHostBufPtr + WRITE_BUF_SIZE;
	    seqWriteInfo->m_hostBufferFilePos[2] = seqWriteInfo->m_hostBufferFilePos[2] + WRITE_BUF_SIZE;
	    break;
	default:
	    velvetLog("Unknown CnySeq host buffer state %d\n", seqWriteInfo->m_hostBuffersInUse);
	    exit(1);
	    break;
    }
}

static void moveCnySeqNucleotides(SequencesWriter *seqWriteInfo)
{
    // move nucleotides in buffer to allow a four byte length value
    // the current sequence may span two buffers

    uint64_t bufIdx = (seqWriteInfo->m_hostBuffersInUse == 2) ? (seqWriteInfo->m_pHostBufPtr - seqWriteInfo->m_pWriteBuffer[1] + WRITE_BUF_SIZE) : (seqWriteInfo->m_pHostBufPtr - seqWriteInfo->m_pWriteBuffer[0]);
    if (bufIdx + 4 > 2 * WRITE_BUF_SIZE) {
	velvetLog("CnySeq bufIdx %ld too large\n", bufIdx);
	exit(1);
    }

    if (bufIdx + 4 >= WRITE_BUF_SIZE) {
	// continue writing to buf[1]
	seqWriteInfo->m_pHostBufPtr = seqWriteInfo->m_pWriteBuffer[1] + (bufIdx + 4 - WRITE_BUF_SIZE);
	seqWriteInfo->m_pHostBufPtrMax = seqWriteInfo->m_pWriteBuffer[1] + WRITE_BUF_SIZE;
	seqWriteInfo->m_hostBufferFilePos[1] = seqWriteInfo->m_hostBufferFilePos[0] + WRITE_BUF_SIZE;
	seqWriteInfo->m_hostBuffersInUse = 2;
    } else
	seqWriteInfo->m_pHostBufPtr += 4;

    seqWriteInfo->m_insertCurrentIndex += 16;

    uint64_t cnt;
    for (cnt = (seqWriteInfo->m_insertLength+3)>>2; cnt > 0; cnt -= 1) {
	seqWriteInfo->m_pWriteBuffer[(bufIdx+4) >> WRITE_BUF_SHFT][(bufIdx+4) & WRITE_BUF_MASK] = seqWriteInfo->m_pWriteBuffer[bufIdx >> WRITE_BUF_SHFT][bufIdx & WRITE_BUF_MASK];
	bufIdx -= 1;
    }
}

static void writeCnySeqNucleotide(uint8_t nucleotide, SequencesWriter *seqWriteInfo)
{
    if (seqWriteInfo->m_insertLength == SHORT_NUCL_LENGTH-1) {
	moveCnySeqNucleotides(seqWriteInfo);
    }
    if ((seqWriteInfo->m_insertCurrentIndex & 0x3) == 0)
	*seqWriteInfo->m_pHostBufPtr = 0;

    *seqWriteInfo->m_pHostBufPtr = *seqWriteInfo->m_pHostBufPtr | (nucleotide << ((seqWriteInfo->m_insertCurrentIndex & 0x3) * 2));

    seqWriteInfo->m_insertLength += 1;
    seqWriteInfo->m_insertCurrentIndex += 1;

    if ((seqWriteInfo->m_insertCurrentIndex & 0x3) == 0) {
	seqWriteInfo->m_pHostBufPtr += 1;

	if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax)
	    cnySeqHostBufferFull(seqWriteInfo);
    }
}

void cnySeqInsertNucleotideString(const char *pReadBuf, SequencesWriter *seqWriteInfo) {
    uint8_t		nucleotide;

    static boolean bInit = false;
    static uint8_t charMap[256];

    if (!bInit) {
	bInit = true;
	// anything unusual defaults to A
	memset(charMap, ADENINE, 256);
	charMap[(int)'C'] = charMap[(int)'c'] = CYTOSINE;
	charMap[(int)'G'] = charMap[(int)'g'] = GUANINE;
	charMap[(int)'T'] = charMap[(int)'t'] = THYMINE;
	charMap[(int)'\0'] = 4;
    }

    for (;;) {
	nucleotide = charMap[(int)*pReadBuf];
	if (nucleotide < 4) {
	    writeCnySeqNucleotide(nucleotide, seqWriteInfo);
	    pReadBuf += 1;
	    continue;
	} else if (nucleotide == 4) {
	    return;
	} else {
	    velvetLog("CnySeq unexpected char %c (%d)\n", *pReadBuf, (int) *pReadBuf);
	    exit(1);
	}
    }
}

SequencesWriter * openCnySeqForWrite(const char *unifiedSeqFileName)
{
    SequencesWriter *seqWriteInfo = callocOrExit(1, SequencesWriter);
    seqWriteInfo->m_pWriteBuffer[0] = NULL;
    seqWriteInfo->m_pWriteBuffer[1] = NULL;
    seqWriteInfo->m_pWriteBuffer[2] = NULL;
    char seqNamesFileName[5000];

    strcpy(seqNamesFileName, unifiedSeqFileName);
    strcat(seqNamesFileName, ".names");

#ifdef COLOR
    seqWriteInfo->m_unifiedSeqFileHeader.m_bColor = true;
#else
    seqWriteInfo->m_unifiedSeqFileHeader.m_bColor = false;
#endif

    if ((seqWriteInfo->m_pFile = fopen(unifiedSeqFileName, "wb")) == 0) {
	velvetLog("Unable to open %s for writing\n", unifiedSeqFileName);
	exit(1);
    }

    if ((seqWriteInfo->m_nameFile = fopen(seqNamesFileName, "w")) == 0) {
        velvetLog("Unable to open %s for writing\n", seqNamesFileName);
        exit(1);
    }

    memcpy(&seqWriteInfo->m_unifiedSeqFileHeader.m_magic, "CSQ0", 4);
    seqWriteInfo->m_unifiedSeqFileHeader.m_timeStamp = time(0);
    seqWriteInfo->m_unifiedSeqFileHeader.m_bFileWriteCompleted = false;

    if (fwrite(&seqWriteInfo->m_unifiedSeqFileHeader, sizeof(CnyUnifiedSeqFileHeader), 1, seqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write file %s\n", unifiedSeqFileName);
	exit(1);
    }

    seqWriteInfo->m_insertCurrentIndex = 0;
    seqWriteInfo->m_pWriteBuffer[0] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);
    seqWriteInfo->m_pWriteBuffer[1] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);
    seqWriteInfo->m_pWriteBuffer[2] = mallocOrExit(WRITE_BUF_SIZE, uint8_t);

    seqWriteInfo->m_hostBufferFilePos[0] = sizeof(CnyUnifiedSeqFileHeader);

    seqWriteInfo->m_pHostBufPtr = seqWriteInfo->m_pWriteBuffer[0];
    seqWriteInfo->m_pHostBufPtrMax = seqWriteInfo->m_pWriteBuffer[0] + WRITE_BUF_SIZE;
    seqWriteInfo->m_hostBuffersInUse = 1;
    seqWriteInfo->m_fileSegmentWriteIdx = 0;	// file segment currently being written
    seqWriteInfo->m_unifiedSeqFileHeader.m_sequenceCnt = 0;
    seqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen = ~0LL;
    seqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen = 0;
    seqWriteInfo->m_unifiedSeqFileHeader.m_totalSeqLen = 0;
    return seqWriteInfo;
}

static void alignCnySeqToNextByteBoundary(SequencesWriter *seqWriteInfo)
{
    if ((seqWriteInfo->m_insertCurrentIndex & 0x3) != 0)
	seqWriteInfo->m_pHostBufPtr += 1;

    seqWriteInfo->m_insertCurrentIndex = (seqWriteInfo->m_insertCurrentIndex + 3) & ~0x3LL;

    if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax) {
	cnySeqHostBufferFull(seqWriteInfo);
    }
}

static void writeCnySeqUint8(uint8_t uint8, SequencesWriter *seqWriteInfo)
{
    *seqWriteInfo->m_pHostBufPtr++ = uint8;
    seqWriteInfo->m_insertCurrentIndex += 4;

    if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax) {
	cnySeqHostBufferFull(seqWriteInfo);
    }
}

static void writeCnySeqUint32(uint32_t uint32, SequencesWriter *seqWriteInfo)
{
    int i;
    for (i = 0; i < 4; i += 1)
	writeCnySeqUint8((uint8_t)(uint32 >> (i*8)), seqWriteInfo);
}

void inputCnySeqFileStart(Category category, SequencesWriter *seqWriteInfo)
{
    if (category > REFERENCE) {
	velvetLog("Found category %d beyond max of %d\n", category, REFERENCE);
	exit(1);
    }

    alignCnySeqToNextByteBoundary(seqWriteInfo);
    writeCnySeqUint8(0xc0, seqWriteInfo);
    writeCnySeqUint32(category, seqWriteInfo);
}

void cnySeqInsertStart(SequencesWriter *seqWriteInfo)
{
    seqWriteInfo->m_unifiedSeqFileHeader.m_sequenceCnt += 1;

    alignCnySeqToNextByteBoundary(seqWriteInfo);

    seqWriteInfo->m_insertLength = 0;
    seqWriteInfo->m_pHostLengthBufPtr = seqWriteInfo->m_pHostBufPtr;
    seqWriteInfo->m_pHostLengthBufPtrMax = seqWriteInfo->m_pHostBufPtrMax;
    seqWriteInfo->m_pHostBufPtr += 1;
    seqWriteInfo->m_insertLengthIndex = seqWriteInfo->m_insertCurrentIndex >> 2;	// byte index
    seqWriteInfo->m_insertCurrentIndex += 4; // allow for single byte header
    seqWriteInfo->m_insertStartIndex = seqWriteInfo->m_insertCurrentIndex;

    if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax)
    {
	cnySeqHostBufferFull(seqWriteInfo);
    }

    seqWriteInfo->m_position = 0;
    seqWriteInfo->m_openMask = false;
}

void cnySeqInsertSequenceName(const char *name, IDnum readID, SequencesWriter *seqWriteInfo, Category cat) {
	if (fprintf(seqWriteInfo->m_nameFile, "%s\t%li\t%li\n", name, (long) readID, (long) cat) < 0) {
		velvetLog("Unable to write in name file\n");
		exit(1);
	}
}

void cnySeqInsertReferenceMask(SequencesWriter *seqWriteInfo, Mask *referenceMask) {
	Mask *tmp;
	for (tmp = referenceMask; tmp; tmp = tmp->next) {
		if (fprintf(seqWriteInfo->m_nameFile, "%li\t%li\n", (long) tmp->start, (long) tmp->finish) < 0) {
			velvetLog("Unable to write ref in name file\n");
			exit(1);
		}
	}
}

void cnySeqInsertEnd(SequencesWriter *seqWriteInfo)
{
    uint8_t *tmp;

    // fill last few empty nucleotides with a fixed pattern for consistency checking
    if ((seqWriteInfo->m_insertCurrentIndex & 0x3) != 0) {
	*seqWriteInfo->m_pHostBufPtr |= 0xAA << ((seqWriteInfo->m_insertCurrentIndex & 0x3)*2);
    }

    // collect read length statistics
    if (seqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen > seqWriteInfo->m_insertLength)
	seqWriteInfo->m_unifiedSeqFileHeader.m_minSeqLen = seqWriteInfo->m_insertLength;
    if (seqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen < seqWriteInfo->m_insertLength)
	seqWriteInfo->m_unifiedSeqFileHeader.m_maxSeqLen = seqWriteInfo->m_insertLength;

    seqWriteInfo->m_unifiedSeqFileHeader.m_totalSeqLen += seqWriteInfo->m_insertLength;

    if (seqWriteInfo->m_insertLength >= SHORT_NUCL_LENGTH || seqWriteInfo->m_bIsRef) {
	if (seqWriteInfo->m_bIsRef) {
	    alignCnySeqToNextByteBoundary(seqWriteInfo);

	    if (seqWriteInfo->m_insertLength < SHORT_NUCL_LENGTH) {
		    // the align above points to next byte,
		    // need to back up to last byte of nucl seq
		    seqWriteInfo->m_pHostBufPtr -= 1;
		    moveCnySeqNucleotides(seqWriteInfo);
		    // move to next byte
		    seqWriteInfo->m_pHostBufPtr += 1;
	    }

	    // write out map info
	    int idx;
	    for (idx = 0; idx < 4; idx += 1) {
		if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax)
		    cnySeqHostBufferFull(seqWriteInfo);
		*seqWriteInfo->m_pHostBufPtr++ = (seqWriteInfo->m_refCnt >> (idx*8)) & 0xff;
		seqWriteInfo->m_insertCurrentIndex += 4; // single byte
	    }
	    int refIdx;
	    RefInfoList *refElem = seqWriteInfo->m_refInfoHead;
	    RefInfoList *prev = NULL;
	    for (refIdx = 0; refIdx < seqWriteInfo->m_refCnt; refIdx++) {

		    if (refElem == NULL) {
			    velvetLog("reference but element %d NULL\n", refIdx);
			    exit(1);
		    }
		    alignCnySeqToNextByteBoundary(seqWriteInfo);

		    for (idx = 0; idx < 4; idx += 1) {
			    if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax)
				    cnySeqHostBufferFull(seqWriteInfo);
			    *seqWriteInfo->m_pHostBufPtr++ = (refElem->m_elem.m_referenceID >> (idx*8)) & 0xff;
			    seqWriteInfo->m_insertCurrentIndex += 4; // single byte
		    }
		    for (idx = 0; idx < 4; idx += 1) {
			    if (seqWriteInfo->m_pHostBufPtr == seqWriteInfo->m_pHostBufPtrMax)
				    cnySeqHostBufferFull(seqWriteInfo);
			    *seqWriteInfo->m_pHostBufPtr++ = (refElem->m_elem.m_pos >> (idx*8)) & 0xff;
			    seqWriteInfo->m_insertCurrentIndex += 4; // single byte
		    }

		    prev = refElem;
		    refElem = refElem->next;
		    free(prev);
	    }
	    if (refElem != NULL) {
		    velvetLog("more than %d elements in ref\n", seqWriteInfo->m_refCnt);
		    exit(1);
	    }
	    seqWriteInfo->m_bIsRef = false;
	    seqWriteInfo->m_refInfoHead = NULL;
	    seqWriteInfo->m_refCnt = 0;	    // set ref bit
	    *(seqWriteInfo->m_pHostLengthBufPtr++) = 0xa0 | ((seqWriteInfo->m_insertLength >> 32) & 0x1f);
	} else {
	    // one byte control, four byte length
	    *(seqWriteInfo->m_pHostLengthBufPtr++) = 0x80 | ((seqWriteInfo->m_insertLength >> 32) & 0x1f);
	}
	int idx;
	for (idx = 0; idx < 4; idx += 1) {
	    if (seqWriteInfo->m_pHostLengthBufPtr == seqWriteInfo->m_pHostLengthBufPtrMax) {
		if (seqWriteInfo->m_hostBuffersInUse < 2) {
		    velvetLog("CnySeq m_hostBuffersInUse %d\n", seqWriteInfo->m_hostBuffersInUse);
		    exit(1);
		}
		seqWriteInfo->m_pHostLengthBufPtr = seqWriteInfo->m_pWriteBuffer[1];
	    }
	    *seqWriteInfo->m_pHostLengthBufPtr++ = (seqWriteInfo->m_insertLength >> (idx*8)) & 0xff;
	}

    } else {
	// one byte length;
	*seqWriteInfo->m_pHostLengthBufPtr = (uint8_t)seqWriteInfo->m_insertLength;
    }

    alignCnySeqToNextByteBoundary(seqWriteInfo);

    switch (seqWriteInfo->m_hostBuffersInUse) {
	case 1:	// buf[0] is being written
	    break;
	case 2: // buf[0] and buf[1] are being written, write buf[0] to disk
	    if (fseek(seqWriteInfo->m_pFile, seqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }

	    if (fwrite(seqWriteInfo->m_pWriteBuffer[0], WRITE_BUF_SIZE, 1, seqWriteInfo->m_pFile) != 1) {
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    // swap buf[0] and buf[1]
	    tmp = seqWriteInfo->m_pWriteBuffer[0];
	    seqWriteInfo->m_pWriteBuffer[0] = seqWriteInfo->m_pWriteBuffer[1];
	    seqWriteInfo->m_pWriteBuffer[1] = tmp;

	    seqWriteInfo->m_hostBufferFilePos[0] = seqWriteInfo->m_hostBufferFilePos[1];

	    seqWriteInfo->m_hostBuffersInUse = 1;
	    break;
	case 3: // buf[0], [1] and [2] are in use, write buf[0] and [1] to disk
	    if (fseek(seqWriteInfo->m_pFile, seqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fwrite(seqWriteInfo->m_pWriteBuffer[0], WRITE_BUF_SIZE, 1, seqWriteInfo->m_pFile) != 1){
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fseek(seqWriteInfo->m_pFile, seqWriteInfo->m_hostBufferFilePos[1], SEEK_SET) < 0) {
		velvetLog("Unable to seek in CnyUnifiedSeq\n");
		exit(1);
	    }


	    if (fwrite(seqWriteInfo->m_pWriteBuffer[1], WRITE_BUF_SIZE, 1, seqWriteInfo->m_pFile) != 1){
		velvetLog("Unable to write CnyUnifiedSeq\n");
		exit(1);
	    }

	    // swap buf[0] and buf[2]
	    tmp = seqWriteInfo->m_pWriteBuffer[0];
	    seqWriteInfo->m_pWriteBuffer[0] = seqWriteInfo->m_pWriteBuffer[2];
	    seqWriteInfo->m_pWriteBuffer[2] = tmp;

	    seqWriteInfo->m_hostBufferFilePos[0] = seqWriteInfo->m_hostBufferFilePos[2];

	    seqWriteInfo->m_hostBuffersInUse = 1;
	    break;
    }

    // if ref masks, write mapping info to the names file
    if (seqWriteInfo->m_referenceMask && *(seqWriteInfo->m_referenceMask)) {
	    cnySeqInsertReferenceMask(seqWriteInfo, *(seqWriteInfo->m_referenceMask));
	    // free memory and clear list
	    if (seqWriteInfo->m_maskMemory) {
		    destroyRecycleBin(seqWriteInfo->m_maskMemory);
		    seqWriteInfo->m_maskMemory = NULL;
	    }
	    *(seqWriteInfo->m_referenceMask) = NULL;
    }
}

void closeCnySeqForWrite(SequencesWriter *seqWriteInfo)
{
    // should be only one buffer in use
    if (seqWriteInfo->m_hostBuffersInUse != 1) {
	velvetLog("CnySeq host buffers in use %d\n", seqWriteInfo->m_hostBuffersInUse);
	exit(1);
    }

    if (fseek(seqWriteInfo->m_pFile, seqWriteInfo->m_hostBufferFilePos[0], SEEK_SET) < 0) {
	velvetLog("Unable to seek CnySeq\n");
	exit(1);
    }

    if (fwrite(seqWriteInfo->m_pWriteBuffer[0], (uint32_t)(seqWriteInfo->m_pHostBufPtr - seqWriteInfo->m_pWriteBuffer[0]), 1, seqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write CnySeq\n");
	exit(1);
    }

    seqWriteInfo->m_unifiedSeqFileHeader.m_bFileWriteCompleted = true;
    seqWriteInfo->m_unifiedSeqFileHeader.m_seqNuclStoreSize = seqWriteInfo->m_insertCurrentIndex >> 2;
    seqWriteInfo->m_unifiedSeqFileHeader.m_numCategories = CATEGORIES;

    if (fseek(seqWriteInfo->m_pFile, 0, SEEK_SET) < 0) {
	velvetLog("Unable to seek CnySeq\n");
	exit(1);
    }

    if (fwrite(&seqWriteInfo->m_unifiedSeqFileHeader, sizeof(CnyUnifiedSeqFileHeader), 1, seqWriteInfo->m_pFile) != 1) {
	velvetLog("Unable to write CnySeq\n");
	exit(1);
    }

    if (fclose(seqWriteInfo->m_pFile) < 0) {
	velvetLog("Unable to close CnySeq\n");
	exit(1);
    }

    if (fclose(seqWriteInfo->m_nameFile) < 0) {
	velvetLog("Unable to close names file\n");
	exit(1);
    }

    if (seqWriteInfo->m_pWriteBuffer[0])
        free(seqWriteInfo->m_pWriteBuffer[0]);
    if (seqWriteInfo->m_pWriteBuffer[1])
        free(seqWriteInfo->m_pWriteBuffer[1]);
    if (seqWriteInfo->m_pWriteBuffer[2])
        free(seqWriteInfo->m_pWriteBuffer[2]);
}

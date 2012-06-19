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
#ifndef _BINARYSEQUENCES_H_
#define _BINARYSEQUENCES_H_

#include "recycleBin.h"

struct refInfo_st {
	int32_t         m_referenceID;
	int32_t         m_pos;
};
typedef struct refInfo_st RefInfo;


typedef struct {
	Category	m_numCategories;
	uint32_t	m_magic;
	boolean		m_bDoubleStrand;
	boolean		m_bColor;
	uint64_t	m_sequenceCnt;
	uint64_t	m_timeStamp;
	uint64_t	m_seqNuclStoreSize;
	uint64_t	m_minSeqLen;
	uint64_t	m_maxSeqLen;
	uint64_t	m_totalSeqLen;
	boolean		m_bFileWriteCompleted;
} CnyUnifiedSeqFileHeader;

// Reading
struct sequencesReader_st {
	char *		m_seqFilename;
        FILE *		m_pFile;
	ReadSet *       m_sequences;
	boolean		m_bIsBinary;
	// following are only used by the binary code
	char *		m_namesFilename;
        CnyUnifiedSeqFileHeader m_unifiedSeqFileHeader;
	Category	m_numCategories;
	Category	m_currCategory;
	uint64_t	m_minSeqLen;
	uint64_t	m_maxSeqLen;
	uint8_t *	m_pReadBuffer;
	uint8_t *	m_pReadBufEnd;
	uint64_t	m_readBufPos;
	uint8_t *	m_pCurrentReadPtr;
	uint8_t *	m_pNextReadPtr;
	uint64_t	m_currentNuclReadIdx;
	uint64_t	m_currentReadLength;
	uint32_t        m_refCnt;
	boolean         m_bIsRef;
	char **		m_ppCurrString;
};

#define USF_READ_BUF_SIZE				(64*1024)

ReadSet *importCnyReadSet(char *filename);
void getCnySeqNucl(SequencesReader *seqReadInfo, uint8_t *sequence);
uint32_t readCnySeqUint32(SequencesReader *seqReadInfo);
boolean advanceCnySeqCurrentRead(SequencesReader *seqReadInfo);
FILE *openCnySeqForRead(const char *fileName, CnyUnifiedSeqFileHeader *seqFileHeader);
void resetCnySeqCurrentRead(SequencesReader *seqReadInfo);

struct refInfoList_st {
	RefInfo		m_elem;
	struct refInfoList_st *next;
};
typedef struct refInfoList_st RefInfoList;

// Writing
struct sequencesWriter_st {
        FILE *		m_pFile;
	// following are only used by the binary code
	FILE *          m_nameFile;
        CnyUnifiedSeqFileHeader m_unifiedSeqFileHeader;
	uint64_t	m_insertStartIndex;
	uint64_t	m_insertLength;
	uint64_t	m_insertLengthIndex;
	uint64_t	m_insertCurrentIndex;
	uint32_t	m_hostBuffersInUse;
	uint32_t	m_fileSegmentWriteIdx;
	uint8_t	*	m_pWriteBuffer[3];
	uint8_t *	m_pHostBufPtr;
	uint8_t *	m_pHostLengthBufPtr;
	uint8_t *	m_pHostLengthBufPtrMax;
	uint8_t *	m_pHostBufPtrMax;
	int64_t		m_hostBufferFilePos[3];
	RefInfoList *	m_refInfoHead;
	uint32_t        m_refCnt;
	boolean         m_bIsRef;
	Mask **         m_referenceMask;
	Mask *          m_current;
	Coordinate      m_position;
	boolean         m_openMask;
	RecycleBin *    m_maskMemory;
};
typedef struct sequencesWriter_st SequencesWriter;

SequencesWriter * openCnySeqForWrite(const char *unifiedSeqFileName);
void cnySeqInsertSequenceName(const char *name, IDnum readID, SequencesWriter *seqWriteInfo, Category cat);
void cnySeqInsertReferenceMask(SequencesWriter *seqWriteInfo, Mask *referenceMask);
void cnySeqInsertNucleotideString(const char *pReadBuf, SequencesWriter *seqWriteInfo);
void inputCnySeqFileStart(Category category, SequencesWriter *seqWriteInfo);
void cnySeqInsertStart(SequencesWriter *seqWriteInfo);
void cnySeqInsertEnd(SequencesWriter *seqWriteInfo);
void closeCnySeqForWrite(SequencesWriter *seqWriteInfo);

#endif

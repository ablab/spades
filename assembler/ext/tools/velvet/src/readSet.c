/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

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
#include <ctype.h>

#include "globals.h"
#include "tightString.h"
#include "readSet.h"
#include "utility.h"
#include "binarySequences.h"
#include "autoOpen.h"
#include "kseq.h"

#if !defined(BUNDLEDZLIB)
#include <zlib.h>
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
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

static Mask *allocateMask(SequencesWriter *seqWriteInfo)
{
	if (seqWriteInfo->m_maskMemory == NULL)
		seqWriteInfo->m_maskMemory = newRecycleBin(sizeof(Mask), 10000);

	return (Mask *) allocatePointer(seqWriteInfo->m_maskMemory);
}

static Mask * newMask(SequencesWriter *seqWriteInfo, Coordinate position)
{
	Mask * mask = allocateMask(seqWriteInfo);
	mask->start = position;
	mask->finish = position;
	mask->next = NULL;
	return mask;
}

//
// cmd line args can override the createBinary flag
// note that createBinary is only used by velveth
//
boolean createBinary = false;
inline boolean isCreateBinary()
{
	return createBinary;
}

void setCreateBinary(boolean val)
{
	createBinary = val;
}

ReadSet *newReadSet()
{
	ReadSet *rs = callocOrExit(1, ReadSet);
	return rs;
}

//////////////////////////////////////////////////////////////////////////
// Reference identifiers
//////////////////////////////////////////////////////////////////////////

typedef struct referenceCoordinate_st ReferenceCoordinate;
static Coordinate reference_coordinate_double_strand = true;

struct referenceCoordinate_st {
	char * name;
	Coordinate start;
	Coordinate finish;
	IDnum referenceID;
	IDnum counter;
	boolean positive_strand;
}  ATTRIBUTE_PACKED;

static int compareRefCoords(const void * ptrA, const void * ptrB) {
	ReferenceCoordinate * A = (ReferenceCoordinate *) ptrA;
	ReferenceCoordinate * B = (ReferenceCoordinate *) ptrB;
	int comp = strcmp(A->name, B->name);

	if (comp != 0)
		return comp;
	else if (!reference_coordinate_double_strand && A->positive_strand != B->positive_strand)
		return A->positive_strand > B->positive_strand;
	else {
		if (A->finish > -1 && A->finish < B->start)
			return -1;
		else if (B->finish > -1 && A->start > B->finish) 
			return 1;
		else return 0;
	}
}

typedef struct referenceCoordinateTable_st ReferenceCoordinateTable;

struct referenceCoordinateTable_st {
	ReferenceCoordinate * array;
	IDnum arrayLength;
}  ATTRIBUTE_PACKED;

static ReferenceCoordinateTable * newReferenceCoordinateTable() {
	ReferenceCoordinateTable * table = callocOrExit(1, ReferenceCoordinateTable);
	table->array = NULL;
	table->arrayLength = 0;
	return table;
}

static void printReferenceCoordinateTableStats(ReferenceCoordinateTable * table) {
	IDnum index;
	IDnum counter = 0;

	velvetLog("Reference mapping counters\n");
	velvetLog("Name\tRead mappings\n");

	for (index = 0; index < table->arrayLength; index++) {
		velvetLog("%s\t%li\n", table->array[index].name, (long) table->array[index].counter);
		counter += table->array[index].counter;
	}

	if (counter == 0) {
		velvetLog("WARNING: None of your read mappings recognized the reference sequence!\n");
		velvetLog("Double check that the names are identical between reference fasta headers and SAM/BAM sequences.\n");
	}
}

static void destroyReferenceCoordinateTable(ReferenceCoordinateTable * table) {
	IDnum index;

	if (table->array) {
		printReferenceCoordinateTableStats(table);
		for (index = 0; index < table->arrayLength; index++)
			free(table->array[index].name);
		free(table->array);
	}
	free(table);
}

static void resizeReferenceCoordinateTable(ReferenceCoordinateTable * table, IDnum extraLength) {
	if (table->array == NULL)
		table->array = callocOrExit(extraLength, ReferenceCoordinate);
	else 
		table->array = reallocOrExit(table->array, table->arrayLength + extraLength, ReferenceCoordinate);
}

static ReferenceCoordinate * findReferenceCoordinate(ReferenceCoordinateTable * table, char * name, Coordinate start, Coordinate finish, boolean positive_strand) {
	ReferenceCoordinate * array = table->array;
	ReferenceCoordinate refCoord;
	Coordinate leftIndex = 0;
	Coordinate rightIndex = table->arrayLength - 1;
	Coordinate middleIndex;

	refCoord.name = name;
	refCoord.start = start;
	refCoord.finish = finish;
	refCoord.referenceID = 0;
	refCoord.positive_strand = positive_strand;

	while (true) {
		middleIndex = (rightIndex + leftIndex) / 2;

		if (leftIndex > rightIndex)
			return NULL;
		else if (compareRefCoords(&(array[middleIndex]), &refCoord) == 0)
			return &(array[middleIndex]);
		else if (leftIndex == middleIndex)
			return NULL;
		else if (compareRefCoords(&(array[middleIndex]), &refCoord) > 0)
			rightIndex = middleIndex;
		else
			leftIndex = middleIndex;
	}
}

static void addReferenceCoordinate(ReferenceCoordinateTable * table, char * name, Coordinate start, Coordinate finish, boolean positive_strand) {
	ReferenceCoordinate * refCoord;

	if ((refCoord = findReferenceCoordinate(table, name, start, finish, positive_strand))) {
		velvetLog("Overlapping reference coordinates:\n");
		velvetLog("%s:%lli-%lli\n", name, (long long) start, (long long) finish);
		velvetLog("%s:%lli-%lli\n", refCoord->name, (long long) refCoord->start, (long long) refCoord->finish);
		velvetLog("Exiting...\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	
	refCoord = &(table->array[table->arrayLength++]);

	refCoord->name = name;
	refCoord->start = start;
	refCoord->finish = finish;
	refCoord->referenceID = table->arrayLength;
	refCoord->positive_strand = positive_strand;
	refCoord->counter = 0;
}

static void sortReferenceCoordinateTable(ReferenceCoordinateTable * table) {
	qsort(table->array, table->arrayLength, sizeof(ReferenceCoordinate), compareRefCoords);
}

//////////////////////////////////////////////////////////////////////////
// File reading 
//////////////////////////////////////////////////////////////////////////

static void velvetifySequence(char * str, SequencesWriter *seqWriteInfo) {
	int i;
	char c;
	size_t length = strlen(str);

	for (i = 0; i < length; i++) {
		c = str[i];
		switch (c) {
		case '\n':
		case '\r':
		case EOF:
			str[i] = '\0';
			break;
		case 'A':
		case 'a':
			str[i] = 'A';
			break;
		case 'C':
		case 'c':
			str[i] = 'C';
			break;
		case 'G':
		case 'g':
			str[i] = 'G';
			break;
		case 'T':
		case 't':
			str[i] = 'T';
			break;
		default:
			str[i] = 'N';
		}
		// non NULL indicates ref masks are being created
		if (seqWriteInfo->m_referenceMask != NULL) {
			if (str[i] == 'N') {
				if (seqWriteInfo->m_openMask) {
					seqWriteInfo->m_current->finish++;
				} else if (*(seqWriteInfo->m_referenceMask) == NULL) {
					*(seqWriteInfo->m_referenceMask) = newMask(seqWriteInfo, seqWriteInfo->m_position);
					seqWriteInfo->m_current = *(seqWriteInfo->m_referenceMask);
				} else {
					seqWriteInfo->m_current->next = newMask(seqWriteInfo, seqWriteInfo->m_position);
					seqWriteInfo->m_current = seqWriteInfo->m_current->next;		
				}
				seqWriteInfo->m_openMask = true;
				seqWriteInfo->m_position += 1;
			} else if (str[i] != '\0') {
				seqWriteInfo->m_openMask = false;
				seqWriteInfo->m_position += 1;
			}
		}
	} 
}

static void reverseComplementSequence(char * str)
{
	size_t length = strlen(str);
	size_t i;

	for (i = 0; i < length-1 - i; i++) {
		char c = str[i];
		str[i] = str[length-1 - i];
		str[length-1 - i] = c;
	}

#ifndef COLOR
	for (i = 0; i < length; i++) {
		switch (str[i]) {
		case 'A':
		case 'a':
			str[i] = 'T';
			break;
		case 'C':
		case 'c':
			str[i] = 'G';
			break;
		case 'G':
		case 'g':
			str[i] = 'C';
			break;
		// As in velvetifySequence(), anything unusual ends up as 'A'
		default:
			str[i] = 'A';
			break;
		}
	}
#endif
}

static void writeFastaSequence(FILE * outfile, const char * str)
{
	size_t length = strlen(str);
	size_t start;
	for (start = 0; start < length; start += 60)
		velvetFprintf(outfile, "%.60s\n", &str[start]);
}

void convertSequences(ReadSet * rs)
{
	rs->tSequences = newTightStringArrayFromStringArray(rs->sequences,
							    rs->readCount,
							    &rs->tSeqMem);
	rs->sequences = NULL;
}

// Returns the value of a 32-bit little-endian-stored integer.
static int int32(const unsigned char * ptr)
{
	int x = ptr[3];
	x = (x << 8) | ptr[2];
	x = (x << 8) | ptr[1];
	x = (x << 8) | ptr[0];
	return x;
}

void goToEndOfLine(char *line, FILE * file)
{
	size_t length = strlen(line);
	char c = line[length - 1];

	while (c != '\n')
		c = fgetc(file);
}

static void writeSeqName(char*seq_name, SequencesWriter *seqWriteInfo, Category cat, IDnum *sequenceIndex)
{
	char name[5001];
	if (isCreateBinary()) {
	  	cnySeqInsertStart(seqWriteInfo);
		sprintf(name, ">%s", seq_name);
		cnySeqInsertSequenceName(name, (long) ((*sequenceIndex)++), seqWriteInfo, cat);
	} else {
		velvetFprintf(seqWriteInfo->m_pFile,">%s\t%ld\t%d\n", seq_name, (long) ((*sequenceIndex)++), (int) cat);
	}
}

static void writeSequence(char*seq, SequencesWriter *seqWriteInfo)
{
	char str[100];
	velvetifySequence(seq, seqWriteInfo);
	if (isCreateBinary()) {
		cnySeqInsertNucleotideString(seq, seqWriteInfo);
		cnySeqInsertEnd(seqWriteInfo);
	} else {
		Coordinate start = 0;
		while (start <= strlen(seq)) {
			strncpy(str, seq + start, 60);
			str[60] = '\0';
				velvetFprintf(seqWriteInfo->m_pFile, "%s\n", str);
			start += 60;
		}
	}
}

static void initFastX(SequencesWriter *seqWriteInfo, Category cat)
{
	seqWriteInfo->m_referenceMask = NULL;
	seqWriteInfo->m_position = 0;
	seqWriteInfo->m_openMask = false;

	// Binary file stuff
	if (isCreateBinary() && (cat == REFERENCE)) {
		seqWriteInfo->m_referenceMask = callocOrExit(1, Mask*);
	}
	if (isCreateBinary()) {
		inputCnySeqFileStart(cat, seqWriteInfo);
	}
}

static void cleanupFastX(SequencesWriter *seqWriteInfo, Category cat)
{
	if (seqWriteInfo->m_referenceMask) {
		free(seqWriteInfo->m_referenceMask);
		seqWriteInfo->m_referenceMask = NULL;
	}
}


// Imports sequences from a raw sequence file
// Memory space allocated within this function.
static void readRawFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum * sequenceIndex)
{
	FILE *file;
	const int maxline = 5000;
	char line[5000];
	IDnum counter = 0;

	initFastX(seqWriteInfo, cat);

	if (strcmp(filename, "-"))
		file = fopen(filename, "r");
	else 
		file = stdin;

	if (file != NULL)
		velvetLog("Reading raw file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	while(fgets(line, maxline, file)) { 
		if (strlen(line) >= maxline - 1) {
			velvetLog("Raw sequence files cannot contain reads longer than %i bp\n", maxline - 1);
#ifdef DEBUG
			abort();
#endif
			exit(1);
		}

		writeSeqName("RAW", seqWriteInfo, cat, sequenceIndex);
		writeSequence(line, seqWriteInfo);
		counter++;
	}
	fclose(file);
	cleanupFastX(seqWriteInfo, cat);
	velvetLog("%li reads found.\n", (long) counter);
	velvetLog("Done\n");
}

// Imports sequences from a zipped raw file 
// Memory space allocated within this function.
static void readRawGZFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex)
{
	gzFile file;
	const int maxline = 5000;
	char line[5000];
	IDnum counter = 0;

	initFastX(seqWriteInfo, cat);
	if (strcmp(filename, "-"))
		file = gzopen(filename, "rb");
	else { 
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		velvetLog("Reading zipped raw sequence file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	while(gzgets(file, line, maxline)) {
		if (strlen(line) >= maxline - 1) {
			velvetLog("Raw sequence files cannot contain reads longer than %i bp\n", maxline - 1);
#ifdef DEBUG
			abort();
#endif
			exit(1);
		}

		writeSeqName("RAW", seqWriteInfo, cat, sequenceIndex);
		writeSequence(line, seqWriteInfo);
		counter++;
	}
	gzclose(file);
	cleanupFastX(seqWriteInfo, cat);
	velvetLog("%li reads found.\n", (long) counter);
	velvetLog("Done\n");
}

static void fillReferenceCoordinateTable(char *filename, ReferenceCoordinateTable * refCoords, IDnum counter)
{
	FILE *file;
	const int maxline = 5000;
	char line[5000];
	char * name;
	long long start, finish;
	Coordinate i;
	IDnum index = 0;

	if (strcmp(filename, "-") == 0)
		exitErrorf(EXIT_FAILURE, false, "Cannot read reference sequence from stdin");
	else
		file = fopen(filename, "r");

	if (counter == 0)
		return;

	resizeReferenceCoordinateTable(refCoords,counter);

	while (fgets(line, maxline, file) && index < counter) {
		if (line[0] == '>') {
			name = callocOrExit(strlen(line), char);

			if (strchr(line, ':')) {
				sscanf(strtok(line, ":-\r\n\t "), ">%s", name);
				sscanf(strtok(NULL, ":-\r\n\t "), "%lli", &start);
				sscanf(strtok(NULL, ":-\r\n\t "), "%lli", &finish);
				if (start <= finish)
					addReferenceCoordinate(refCoords, name, start, finish, true);
				else
					addReferenceCoordinate(refCoords, name, finish, start, false);
			} else {
				// Chomping EOL characters and comments
				for (i=strlen(line) - 1; i >= 0; i--)
					if (line[i] == '\n' || line[i] == '\r' || line[i] == ' ' || line[i] == '\t')
						line[i] = '\0';

				strcpy(name, line + 1);
				addReferenceCoordinate(refCoords, name, 1, -1, true);
			}

			index++;
		}
	}

	sortReferenceCoordinateTable(refCoords);
}

#define FASTQ 1
#define FASTA 2
#define FASTA_GZ 5
#define FASTQ_GZ 6
#define SAM 8
#define BAM 9
#define RAW 10
#define RAW_GZ 11
#define AUTO 12

static gzFile openFastXFile(int fileType, char*filename)
{
	gzFile file;
	char c;

    	// Choose file or stdin
	if (strcmp(filename, "-")==0) {
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	} else {
		file = gzopen(filename, "rb");
	}

	// Verify filetype
	c = gzgetc(file);
	switch (fileType) {
	case FASTA:
	case FASTA_GZ:
	  if (c != EOF && c!='>')
	    exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastA format", filename);
	  break;
	case FASTQ:
	case FASTQ_GZ:
	  if (c != EOF && c!='@')
	    exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastQ format", filename);
	  break;
	}
	gzungetc(c, file);


	if (file != NULL) {
		char *type;
		switch (fileType) {
		case FASTA:
		case FASTA_GZ: type = "FastA"; break;
		case FASTQ:
		case FASTQ_GZ: type = "FastQ"; break;
		default: type = ""; break;
		}
		velvetLog("Reading %s file %s;\n", type, filename);
	} else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	return file;
}

typedef struct {
	gzFile gzFile;
	AutoFile *autoFile;
} FileGZOrAuto;

size_t fileGZOrAuto_read(FileGZOrAuto kseq_file, void *ptr, size_t size)
{
	if (kseq_file.gzFile)
        	return gzread(kseq_file.gzFile, ptr, size);
        else
          return fread(ptr, 1, size, kseq_file.autoFile->file);
}

void fileGZOrAuto_close(FileGZOrAuto kseq_file)
{
	if (kseq_file.gzFile)
        	gzclose(kseq_file.gzFile);
        else
          closeFileAuto(kseq_file.autoFile);
}

char const* charToType(char c)
{
	switch(c) {
        case '>': return "FastA";
        case '@': return "FastQ";
        default: return "Unknown";
        }
}

// Define mode to use kseq in
KSEQ_INIT(FileGZOrAuto, fileGZOrAuto_read)

// Read in FastA or FastQ files in compressed or gz format
static void readFastXFile(int fileType, SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum * sequenceIndex, ReferenceCoordinateTable * refCoords)
{
	kseq_t *seq;
	FileGZOrAuto file;
	IDnum counter = 0;

        file.gzFile = file.autoFile = NULL;
        if (fileType == AUTO) {
        	file.autoFile = openFileAuto(filename);
                if (!file.autoFile)
                	exitErrorf(EXIT_FAILURE, false, "Unable to open file '%s' in auto mode", filename);
                velvetLog("Reading file '%s' using '%s' as %s\n", filename, file.autoFile->decompressor, charToType(file.autoFile->first_char));
        } else
          file.gzFile = openFastXFile(fileType, filename);

	initFastX(seqWriteInfo, cat);
	// Read a sequence at a time
	seq = kseq_init(file);
	while (kseq_read(seq) >= 0) {
		counter++;
		writeSeqName(seq->name.s, seqWriteInfo, cat, sequenceIndex);
		writeSequence(seq->seq.s, seqWriteInfo);
	}

	kseq_destroy(seq);
	fileGZOrAuto_close(file);

  	if (cat == REFERENCE) {
		fillReferenceCoordinateTable(filename, refCoords, counter);
	}
	cleanupFastX(seqWriteInfo, cat);

	velvetLog("%li sequences found\n", (long) counter);
	velvetLog("Done\n");
}

static void readFastXPair(int fileType, SequencesWriter *seqWriteInfo, char *filename1, char *filename2, Category cat, IDnum * sequenceIndex)
{
	kseq_t *seq1, *seq2;
	FileGZOrAuto file1, file2;
	IDnum counter = 0;

	if (cat==REFERENCE)
		exitErrorf(EXIT_FAILURE, false, "Cannot read reference sequence in 'separate' read mode");

        file1.gzFile = file1.autoFile = NULL;
        file2.gzFile = file2.autoFile = NULL;
        if (fileType == AUTO) {
        	file1.autoFile = openFileAuto(filename1);
                if (!file1.autoFile)
                	exitErrorf(EXIT_FAILURE, false, "Unable to open file '%s' in auto mode", filename1);
                velvetLog("Reading file '%s' using '%s' as %s\n", filename1, file1.autoFile->decompressor, charToType(file1.autoFile->first_char));
        	file2.autoFile = openFileAuto(filename2);
                if (!file2.autoFile)
                	exitErrorf(EXIT_FAILURE, false, "Unable to open file '%s' in auto mode", filename2);
                velvetLog("Reading file '%s' using '%s' as %s\n", filename2, file2.autoFile->decompressor, charToType(file2.autoFile->first_char));
        } else {
        	file1.gzFile = openFastXFile(fileType, filename1);
                file2.gzFile = openFastXFile(fileType, filename2);
        }
	initFastX(seqWriteInfo, cat);

	// Read a sequence at a time
	seq1 = kseq_init(file1);
	seq2 = kseq_init(file2);
	while (kseq_read(seq1) >= 0) {
		counter++;
		writeSeqName(seq1->name.s, seqWriteInfo, cat, sequenceIndex);
		writeSequence(seq1->seq.s, seqWriteInfo);

		if (kseq_read(seq2) < 0)
		  exitErrorf(EXIT_FAILURE, false, "Right sequence file '%s' has too few sequences", filename2);

		counter++;
		writeSeqName(seq2->name.s, seqWriteInfo, cat, sequenceIndex);
		writeSequence(seq2->seq.s, seqWriteInfo);
	}
	if (kseq_read(seq2) >= 0)
		  exitErrorf(EXIT_FAILURE, false, "Right sequence file '%s' has too many sequences", filename2);

	kseq_destroy(seq1);
	kseq_destroy(seq2);

        fileGZOrAuto_close(file1);
        fileGZOrAuto_close(file2);

	cleanupFastX(seqWriteInfo, cat);

	velvetLog("%li sequences found in total in the paired sequence files\n", (long) counter);
	velvetLog("Done\n");
}

static void addMapping(boolean orientation, Coordinate pos, char * seq, ReferenceCoordinate * refCoord, char * buffer, SequencesWriter * seqWriteInfo, RefInfoList ** refTail, size_t * buffer_size) {
	if (isCreateBinary()) {
		seqWriteInfo->m_bIsRef = true;
		RefInfoList *refElem = callocOrExit(1, RefInfoList);
		if (refCoord->positive_strand) {
			refElem->m_elem.m_referenceID = (long) orientation * refCoord->referenceID;
			refElem->m_elem.m_pos = (long long) (pos - refCoord->start);
		} else {
			refElem->m_elem.m_referenceID = (long) -orientation * refCoord->referenceID;
			refElem->m_elem.m_pos = (long long) (refCoord->finish - pos - strlen(seq));
		}
		refElem->next = NULL;
		if (seqWriteInfo->m_refInfoHead == NULL) {
			seqWriteInfo->m_refInfoHead = refElem;
		} else {
			(*refTail)->next = refElem;
		}
		*refTail = refElem;
		seqWriteInfo->m_refCnt++;
	} else {
		if (refCoord->positive_strand) {
			sprintf(buffer, "%sM\t%li\t%lli\n", buffer, (long) orientation * refCoord->referenceID, (long long) (pos - refCoord->start));
		} else 
			sprintf(buffer, "%sM\t%li\t%lli\n", buffer, (long) - orientation * refCoord->referenceID, (long long) (refCoord->finish - pos - strlen(seq)));

		if (*buffer_size - strlen(buffer) < 100) {
			*buffer_size += 1000;
			buffer = reallocOrExit(buffer, *buffer_size, char);
		}
	}

	// Increment counter
	refCoord->counter++;
}

static void writeMappedSequence(IDnum * sequenceIndex, Category cat, Category prev_cat, char * previous_seq, char * previous_qname, char * previous_qname_pairing, char * buffer, SequencesWriter * seqWriteInfo) {
	char print_qname[5000];
	if (isCreateBinary()) {
		if (prev_cat != cat) {
			inputCnySeqFileStart(cat, seqWriteInfo);
			prev_cat = cat;
		}
		cnySeqInsertStart(seqWriteInfo);
		cnySeqInsertNucleotideString(previous_seq, seqWriteInfo);
		sprintf(print_qname, ">%s%s", previous_qname, previous_qname_pairing);
		cnySeqInsertSequenceName(print_qname, (long) ((*sequenceIndex)++), seqWriteInfo, cat);
		cnySeqInsertEnd(seqWriteInfo);
	} else {
		velvetFprintf(seqWriteInfo->m_pFile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
		(long) ((*sequenceIndex)++), (int) cat);
		writeFastaSequence(seqWriteInfo->m_pFile, previous_seq);
		velvetFprintf(seqWriteInfo->m_pFile, "%s", buffer);
		strcpy(buffer, "");
	}
}

static void readCigar(char * cigar, boolean orientation, Coordinate pos, char * seq, ReferenceCoordinate * refCoord, char * buffer, SequencesWriter * seqWriteInfo, RefInfoList ** refTail, size_t * buffer_size) {
	long long cigar_num;
	int cigar_index;
	char c;

	if (strlen(cigar) == 1 && cigar[0] == '*')
		;
	else {
		cigar_num = 0;
		for (cigar_index = 0; cigar_index < strlen(cigar); cigar_index++) {
			c = cigar[cigar_index];
			if (c == 'M' || c == '=' || c == 'X') {
				if (refCoord->finish < 0 || pos < refCoord->finish)
					addMapping(orientation, pos, seq, refCoord, buffer, seqWriteInfo, refTail, buffer_size);
				cigar_num = 0;
			} else if (c == 'S' || c == 'I') {
				pos -= cigar_num;
				cigar_num = 0;
			} else if (c == 'D' || c == 'N') {
				pos += cigar_num;
				cigar_num = 0;
			} else if (c == 'H' || c == 'P') {
				cigar_num = 0;
			} else if (isdigit(c)) {
				cigar_num = 10 * cigar_num + (c - 48);
			} else {
				abort();
			}
		}
	}
}

static void readSAMFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex, ReferenceCoordinateTable * refCoords)
{
	char line[5000];
	unsigned long lineno;
	IDnum readCount = 0;
	char previous_qname_pairing[10];
	char previous_qname[5000];
	char previous_seq[5000];
	boolean previous_paired = false;
	Category prev_cat = cat;
	Category apparentCat;
	ReferenceCoordinate * refCoord;
	RefInfoList *refTail = NULL;
	seqWriteInfo->m_referenceMask = NULL;   // no ref masks for SAM/BAM
	seqWriteInfo->m_position = 0;
	seqWriteInfo->m_openMask = false;

	size_t buffer_size = 5000;
	char * buffer = callocOrExit(buffer_size, char);

	if (cat == REFERENCE) {
		velvetLog("SAM file %s cannot contain reference sequences.\n", filename);
		velvetLog("Please check the command line.\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}

	FILE *file = (strcmp(filename, "-") != 0)? fopen(filename, "r") : stdin;
	if (file)
		velvetLog("Reading SAM file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);
	if (isCreateBinary()) {
		inputCnySeqFileStart(cat, seqWriteInfo);
	}
	strcpy(previous_qname, "");
	for (lineno = 1; fgets(line, sizeof(line), file); lineno++) {
		if (line[0] != '@') {
			char *qname, *flag, *seq, *rname, *cigar;
			long long pos;
			int orientation;
			int i;

			qname = strtok(line, "\t");
			flag  = strtok(NULL, "\t");
			rname = strtok(NULL, "\t");
			sscanf(strtok(NULL, "\t"), "%lli", &pos);
			orientation = 1;

			// Mapping scor
			(void) strtok(NULL, "\t");
			cigar = strtok(NULL, "\t");

			// Columns 7,8,9 are paired name, position and score
			for (i = 7; i < 10; i++)
				(void) strtok(NULL, "\t");
			seq = strtok(NULL, "\t");

			if (seq == NULL) {
				velvetFprintf(stderr,
					"Line #%lu: ignoring SAM record with too few fields\n",
					lineno);
			}
			else if (strcmp(seq, "*") == 0) {
				velvetFprintf(stderr,
					"Line #%lu: ignoring SAM record with omitted SEQ field\n",
					lineno);
			}
			else {
				// Accept flags represented in either decimal or hex:
				int flagbits = strtol(flag, NULL, 0);

				if (flagbits & 0x4)
				    strcpy(rname, "");

				const char *qname_pairing = "";
				if (flagbits & 0x40)
					qname_pairing = "/1";
				else if (flagbits & 0x80)
					qname_pairing = "/2";

				if (flagbits & 0x10) {
					orientation = -1;
					reverseComplementSequence(seq);
				}

				// Determine if paired to previous read
				boolean same_name = (strcmp(qname, previous_qname) == 0);
				if (readCount && (!same_name || strcmp(qname_pairing, previous_qname_pairing) != 0)) {
					if (cat % 2 && !same_name && !previous_paired)
						apparentCat = cat - 1;
					else
						apparentCat = cat;

					previous_paired = (cat % 2 && same_name);

					writeMappedSequence(sequenceIndex, apparentCat, prev_cat, previous_seq, previous_qname, previous_qname_pairing, buffer, seqWriteInfo);
					prev_cat = apparentCat;
				}

				if (!(flagbits & 0x4) && (refCoord = findReferenceCoordinate(refCoords, rname, (Coordinate) pos, (Coordinate) pos + strlen(seq) - 1, orientation))) {
					readCigar(cigar, orientation, pos, seq, refCoord, buffer, seqWriteInfo, &refTail, &buffer_size);
				}

				strcpy(previous_qname, qname);
				strcpy(previous_qname_pairing, qname_pairing);
				strcpy(previous_seq, seq);
				velvetifySequence(previous_seq, seqWriteInfo);

				readCount++;
			}
		}
	}

	if (readCount) {
		if (cat % 2 && !previous_paired)
			apparentCat = cat - 1;
		else
			apparentCat = cat;
		writeMappedSequence(sequenceIndex, apparentCat, prev_cat, previous_seq, previous_qname, previous_qname_pairing, buffer, seqWriteInfo);
	}

	free(buffer);
	fclose(file);
	velvetLog("%lu reads found.\n", (long) readCount);
	velvetLog("Done\n");
}

static int readBAMint32(gzFile file)
{
	unsigned char buffer[4];
	if (gzread(file, buffer, 4) != 4)
		exitErrorf(EXIT_FAILURE, false, "BAM file header truncated");

	return int32(buffer);
}

static void readBAMFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex, ReferenceCoordinateTable * refCoords)
{
	size_t seqCapacity = 0;
	char *seq = NULL;
	char cigar[5000];
	char cigar_buffer[5000];
	size_t bufferCapacity = 4;
	unsigned char *buffer = mallocOrExit(bufferCapacity, unsigned char);
	unsigned long recno, readCount;
	int i, refCount;
	gzFile file;
	char previous_qname_pairing[10];
	char previous_qname[5000];
	char previous_seq[5000];
	boolean previous_paired = false;
	Category prev_cat = cat;
	Category apparentCat;
	char ** refNames;
	ReferenceCoordinate * refCoord;
	seqWriteInfo->m_referenceMask = NULL;   // no ref masks for SAM/BAM
	seqWriteInfo->m_position = 0;
	seqWriteInfo->m_openMask = false;

	RefInfoList *refTail = NULL;
	size_t mapBuffer_size = 1000;
	char * mapBuffer = callocOrExit(mapBuffer_size, char);

	if (cat == REFERENCE) {
		velvetLog("BAM file %s cannot contain reference sequences.\n", filename);
		velvetLog("Please check the command line.\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}

	if (strcmp(filename, "-") != 0)
		file = gzopen(filename, "rb");
	else {
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		velvetLog("Reading BAM file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	if (! (gzread(file, buffer, 4) == 4 && memcmp(buffer, "BAM\1", 4) == 0))
		exitErrorf(EXIT_FAILURE, false, "%s is not in BAM format", filename);

	// Skip header text
	if (gzseek(file, readBAMint32(file), SEEK_CUR) == -1)
		exitErrorf(EXIT_FAILURE, false, "gzseek failed");

	// Skip header reference list
	refCount = readBAMint32(file);
	refNames = callocOrExit(refCount, char *);
	for (i = 0; i < refCount; i++) {
		int strLength;

		if (gzread(file, buffer, 4) != 4)
			exitErrorf(EXIT_FAILURE, false, "BAM alignment record truncated");

		strLength = int32(buffer);
		refNames[i] = callocOrExit(strLength, char);
		
		if (bufferCapacity < 4 + strLength) {
			bufferCapacity = 4 + strLength + 4096;
			buffer = reallocOrExit(buffer, bufferCapacity, unsigned char);
		}

		if (gzread(file, buffer, 4 + strLength) != 4 + strLength)
			exitErrorf(EXIT_FAILURE, false, "BAM alignment record truncated");

		strcpy(refNames[i], (char *) buffer); 
	}
	if (isCreateBinary()) {
		inputCnySeqFileStart(cat, seqWriteInfo);
	}
	strcpy(previous_qname, "");
	readCount = 0;
	for (recno = 1; gzread(file, buffer, 4) == 4; recno++) {
		int blockSize = int32(buffer);
		int readLength;

		if (bufferCapacity < 4 + blockSize) {
			bufferCapacity = 4 + blockSize + 4096;
			buffer = reallocOrExit(buffer, bufferCapacity, unsigned char);
		}

		if (gzread(file, &buffer[4], blockSize) != blockSize)
			exitErrorf(EXIT_FAILURE, false, "BAM alignment record truncated");

		readLength = int32(&buffer[20]);
		if (readLength == 0) {
			velvetFprintf(stderr,
				"Record #%lu: ignoring BAM record with omitted SEQ field\n",
				recno);
		}
		else {
			int readNameLength = buffer[12];
			int flag_nc = int32(&buffer[16]);
			int flagbits = flag_nc >> 16;
			int cigarLength = flag_nc & 0xffff;
			char *qname = (char *)&buffer[36];
			uint32_t *rawcigar = (uint32_t *) &buffer[36 + readNameLength];
			unsigned char *rawseq =
					&buffer[36 + readNameLength + 4 * cigarLength];
			int rID = int32(&buffer[4]);
			// NOTE: BAM file coords are 0-based, not 1-based like SAM files
			// No comment
			long long pos = int32(&buffer[8]) + 1;
			int orientation = 1;

			const char *qname_pairing = "";
			if (flagbits & 0x40)
				qname_pairing = "/1";
			else if (flagbits & 0x80)
				qname_pairing = "/2";

			strcpy(cigar, "");
			for (i = 0; i < cigarLength; i++) {
				static const char decode_ops[] = "MIDNSHP=X";
				uint32_t packed = *(rawcigar++);
				sprintf(cigar_buffer, "%i%c", packed >> 4, decode_ops[packed & 0xf]);
				strcat(cigar, cigar_buffer);
			}

			if (seqCapacity < readLength + 1) {
				seqCapacity = readLength * 2 + 1;
				seq = reallocOrExit(seq, seqCapacity, char);
			}

			for (i = 0; i < readLength; i += 2) {
				static const char decode_bases[] = "=ACMGRSVTWYHKDBN";
				unsigned int packed = *(rawseq++);
				seq[i] = decode_bases[packed >> 4];
				seq[i+1] = decode_bases[packed & 0xf];
			}
			seq[readLength] = '\0';

			if (flagbits & 0x10) {
				orientation = -1;
				reverseComplementSequence(seq);
			}

			// Determine if paired to previous read
			boolean same_name = (strcmp(qname, previous_qname) == 0);
			if (readCount > 0 && (!same_name || strcmp(qname_pairing, previous_qname_pairing) != 0)) {
				if (cat % 2 && !same_name && !previous_paired)
					apparentCat = cat - 1;
				else
					apparentCat = cat;

				previous_paired = (cat % 2 && same_name);

				writeMappedSequence(sequenceIndex, apparentCat, prev_cat, previous_seq, previous_qname, previous_qname_pairing, mapBuffer, seqWriteInfo);
				prev_cat = apparentCat;
			}

			if (!(flagbits & 0x4) && (refCoord = findReferenceCoordinate(refCoords, refNames[rID], (Coordinate) pos, (Coordinate) pos + strlen(seq) - 1, orientation)))
				readCigar(cigar, orientation, pos, seq, refCoord, mapBuffer, seqWriteInfo, &refTail, &mapBuffer_size);

			strcpy(previous_qname, qname);
			strcpy(previous_qname_pairing, qname_pairing);
			strcpy(previous_seq, seq);
			velvetifySequence(previous_seq, seqWriteInfo);

			readCount++;
		}
	}

	if (readCount) {
		if (cat % 2 && !previous_paired)
			apparentCat = cat - 1;
		else
			apparentCat = cat;
		writeMappedSequence(sequenceIndex, apparentCat, prev_cat, previous_seq, previous_qname, previous_qname_pairing, mapBuffer, seqWriteInfo);
	}

	free(seq);
	free(buffer);
	free(mapBuffer);

	gzclose(file);
	velvetLog("%lu reads found.\n", readCount);
	velvetLog("Done\n");
}


static void printUsage()
{
	puts("Usage:");
	puts("./velveth directory hash_length {[-file_format][-read_type][-separate|-interleaved] filename} [options]");
	puts("");
	puts("\tdirectory\t\t: directory name for output files");
	printf("\thash_length\t\t: odd integer (if even, it will be decremented) <= %i (if above, will be reduced)\n", MAXKMERLENGTH);
	puts("\tfilename\t\t: path to sequence file or - for standard input");	
	puts("");
	puts("File format options:");
	puts("\t-fasta");
	puts("\t-fastq");
	puts("\t-raw");
	puts("\t-fasta.gz");
	puts("\t-fastq.gz");
	puts("\t-raw.gz");
	puts("\t-sam");
	puts("\t-bam");
        puts("\t-fmtAuto");
	puts("");
	puts("Read type options:");
	puts("\t-short");
	puts("\t-shortPaired");
	puts("\t-short2");
	puts("\t-shortPaired2");
	puts("\t-long");
	puts("\t-longPaired");
	puts("\t-reference");
	puts("");
	puts("Options:");
	puts("\t-strand_specific\t: for strand specific transcriptome sequencing data (default: off)");
	puts("");
	puts("Output:");
	puts("\tdirectory/Roadmaps");
	puts("\tdirectory/Sequences");
	puts("\t\t[Both files are picked up by graph, so please leave them there]");
}

// General argument parser for most functions
// Basically a reused portion of toplevel code dumped into here
void parseDataAndReadFiles(char * filename, int argc, char **argv, boolean * double_strand, boolean * noHash)
{
	int argIndex = 1;
	int filetype = FASTA;
	Category cat = 0;
	IDnum sequenceIndex = 1;
	short short_var;
	ReferenceCoordinateTable * refCoords = newReferenceCoordinateTable();
	boolean reuseSequences = false;
	boolean separate_pair_files = false;

	if (argc < 2) {
		printUsage();
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}

	for (argIndex = 1; argIndex < argc; argIndex++) {
		if (strcmp(argv[argIndex], "-strand_specific") == 0) {
			*double_strand = false;
			reference_coordinate_double_strand = false;
		} else if (strcmp(argv[argIndex], "-reuse_Sequences") == 0) {
			reuseSequences = true;
		} else if (strcmp(argv[argIndex], "-noHash") == 0) {
			*noHash = true;
		}
	}

	if (reuseSequences) 
		return;

	SequencesWriter * seqWriteInfo = NULL;
	if (isCreateBinary()) {
		seqWriteInfo = openCnySeqForWrite(filename);
		seqWriteInfo->m_unifiedSeqFileHeader.m_bDoubleStrand = *double_strand;
		// file is already open
	} else {
		seqWriteInfo = callocOrExit(1, SequencesWriter);
		seqWriteInfo->m_pFile = fopen(filename, "w");
	}

	for (argIndex = 1; argIndex < argc; argIndex++) {
		if (argv[argIndex][0] == '-' && strlen(argv[argIndex]) > 1) {

			if (strcmp(argv[argIndex], "-fastq") == 0)
				filetype = FASTQ;
			else if (strcmp(argv[argIndex], "-fasta") == 0)
				filetype = FASTA;
			else if (strcmp(argv[argIndex], "-fastq.gz") == 0)
				filetype = FASTQ_GZ;
			else if (strcmp(argv[argIndex], "-fasta.gz") == 0)
				filetype = FASTA_GZ;
			else if (strcmp(argv[argIndex], "-sam") == 0)
				filetype = SAM;
			else if (strcmp(argv[argIndex], "-bam") == 0)
				filetype = BAM;
			else if (strcmp(argv[argIndex], "-raw") == 0)
				filetype = RAW;
			else if (strcmp(argv[argIndex], "-raw.gz") == 0)
				filetype = RAW_GZ;
			else if (strcmp(argv[argIndex], "-fmtAuto") == 0)
				filetype = AUTO;
			else if (strcmp(argv[argIndex], "-short") == 0)
				cat = 0;
			else if (strcmp(argv[argIndex], "-shortPaired") ==
				 0)
				cat = 1;
			else if (strncmp
				 (argv[argIndex], "-shortPaired",
				  12) == 0) {
				sscanf(argv[argIndex], "-shortPaired%hd", &short_var);
				cat = (Category) short_var;
				if (cat < 1 || cat > CATEGORIES) {
					velvetLog("Unknown option: %s\n",
					       argv[argIndex]);
#ifdef DEBUG 
					abort();
#endif 
					exit(1);
				}
				cat--;
				cat *= 2;
				cat++;
			} else if (strncmp(argv[argIndex], "-short", 6) ==
				   0) {
				sscanf(argv[argIndex], "-short%hd", &short_var);
				cat = (Category) short_var;
				if (cat < 1 || cat > CATEGORIES) {
					velvetLog("Unknown option: %s\n",
					       argv[argIndex]);
#ifdef DEBUG 
					abort();
#endif 
					exit(1);
				}
				cat--;
				cat *= 2;
			} else if (strcmp(argv[argIndex], "-long") == 0)
				cat = LONG;		// CATEGORIES * 2;
			else if (strcmp(argv[argIndex], "-longPaired") == 0)
				cat = LONG_PAIRED;	// CATEGORIES * 2 + 1;
			else if (strcmp(argv[argIndex], "-reference") == 0)
				cat = REFERENCE;	// CATEGORIES * 2 + 2
			else if (strcmp(argv[argIndex], "-strand_specific") == 0) {
				*double_strand = false;
				reference_coordinate_double_strand = false;
			} else if (strcmp(argv[argIndex], "-noHash") == 0) {
				;
			} else if (strcmp(argv[argIndex], "-create_binary") == 0) {
				;
			} else if (strcmp(argv[argIndex], "-interleaved") == 0) {
				separate_pair_files = false;
			} else if (strcmp(argv[argIndex], "-separate") == 0) {
				separate_pair_files = true;
			}
			else {
				velvetLog("Unknown option: %s\n",
				       argv[argIndex]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}

			continue;
		}

		if (cat == -1)
			continue;

		switch (filetype) {
		case FASTA:
		case FASTQ:
		case FASTA_GZ:
		case FASTQ_GZ:
		case AUTO:
			// Separate files for paired reads?  Note odd categories used for paired read type
			if (separate_pair_files && cat%2==1) {
				argIndex++;
				if (argIndex>=argc)
					exitErrorf(EXIT_FAILURE, false, "Require left & right filename for -separate mode");
				readFastXPair(filetype, seqWriteInfo, argv[argIndex-1], argv[argIndex], cat, &sequenceIndex);
			} else {
				readFastXFile(filetype, seqWriteInfo, argv[argIndex], cat, &sequenceIndex, refCoords);
			}
			break;
		case RAW:
			if (separate_pair_files && cat%2==1) {
                        	exitErrorf(EXIT_FAILURE, false, "Currently do not support -separate mode for RAW");
                        }
			readRawFile(seqWriteInfo, argv[argIndex], cat, &sequenceIndex);
			break;
		case RAW_GZ:
			if (separate_pair_files && cat%2==1) {
                        	exitErrorf(EXIT_FAILURE, false, "Currently do not support -separate mode for RAW");
                        }
			readRawGZFile(seqWriteInfo, argv[argIndex], cat, &sequenceIndex);
			break;
		case SAM:
			readSAMFile(seqWriteInfo, argv[argIndex], cat, &sequenceIndex, refCoords);
			break;
		case BAM:
			readBAMFile(seqWriteInfo, argv[argIndex], cat, &sequenceIndex, refCoords);
			break;
		default:
			velvetLog("Screw up in parser... exiting\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(1);
		}
	}

	destroyReferenceCoordinateTable(refCoords);
	if (isCreateBinary()) {
		closeCnySeqForWrite(seqWriteInfo);
	} else {
		fclose(seqWriteInfo->m_pFile);
}
	if (seqWriteInfo) {
	    free(seqWriteInfo);
	}
}

void createReadPairingArray(ReadSet* reads)
{
	IDnum index;
	IDnum *mateReads = mallocOrExit(reads->readCount, IDnum);
	Category cat = 0;
	int phase = 0;

	for (index = 0; index < reads->readCount; index++) 
		mateReads[index] = -1;

	reads->mateReads = mateReads;

	for (index = 0; index < reads->readCount; index++)
	{
		// Paired category
		if (cat & 1)
		{
			// Leaving the paired category
			if (reads->categories[index] != cat)
			{
				if (phase == 1)
				{
					reads->mateReads[index - 1] = -1;
					reads->categories[index - 1]--;
					phase = 0;
				}
				cat = reads->categories[index];
				// Into another paired category
				if (cat & 1)
				{
					reads->mateReads[index] = index + 1;
					phase = 1;
				}
			}
			else if (phase == 0)
			{
				reads->mateReads[index] = index + 1;
				phase = 1;
			}
			else
			{
				reads->mateReads[index] = index - 1;
				phase = 0;
			}
		}
		// Leaving an unpaired category
		else if (reads->categories[index] != cat)
		{
			cat = reads->categories[index];
			// Into a paired category
			if (cat & 1)
			{
				reads->mateReads[index] = index + 1;
				phase = 1;
			}
		}
	}
}

int pairedCategories(ReadSet * reads)
{
	boolean pairedCat[CATEGORIES + 1];
	int pairedCatCount = 0;
	IDnum index;

	for (index = 0; index <= CATEGORIES; index++)
		pairedCat[index] = 0;

	for (index = 0; index < reads->readCount; index++) {
		if (reads->categories[index] & 1 && !pairedCat[reads->categories[index] / 2]) {
			pairedCat[reads->categories[index] / 2] = true;
			if (pairedCatCount++ == CATEGORIES)
				break;
		} 	
	}

	return pairedCatCount;
}

boolean isSecondInPair(ReadSet * reads, IDnum index)
{
	return reads->secondInPair[index / 8] & (1 << (index & 7));
}

void computeSecondInPair(ReadSet * reads)
{
	IDnum index;
	Category currentCat = 0;
	Category previousCat = 0;
	int phase = 0;

	if (reads->secondInPair)
		free (reads->secondInPair);
	reads->secondInPair = callocOrExit((reads->readCount + 7) / 8, unsigned char);

	for (index = 0; index < reads->readCount; index++)
	{
		currentCat = reads->categories[index];
		if (currentCat & 1)
		{
			if (previousCat == currentCat)
			{
				if (phase == 0)
				{
					phase = 1;
				}
				else
				{
					reads->secondInPair[index / 8] |= (1 << (index & 7));
					phase = 0;
				}
			}
			else {
				phase = 1;
				if (index > 0 && previousCat & 1 && !isSecondInPair(reads, index - 1))
				    reads->categories[index - 1] = (reads->categories[index - 1] / 2) * 2;
			}
		}
		previousCat = currentCat;
	}

	// Safeguard against odd sets of reads
	if (!isSecondInPair(reads, reads->readCount - 1)) {
		reads->categories[reads->readCount - 1] = (reads->categories[reads->readCount - 1] / 2) * 2;
	}
}

void detachDubiousReads(ReadSet * reads, boolean * dubiousReads)
{
	IDnum index;
	IDnum pairID;
	IDnum sequenceCount = reads->readCount;
	IDnum *mateReads = reads->mateReads;

	if (dubiousReads == NULL || mateReads == NULL)
		return;

	for (index = 0; index < sequenceCount; index++) {
		if (!dubiousReads[index] || reads->categories[index] % 2 == 0 )
			continue;

		if (isSecondInPair(reads, index))
		    pairID = index - 1;
		else
		    pairID = index + 1;

		reads->categories[index] = (reads->categories[index] / 2) * 2;
		reads->categories[pairID] = (reads->categories[pairID] / 2) * 2;
	}
}

ReadSet *importReadSet(char *filename)
{
	FILE *file = fopen(filename, "r");
	char *sequence = NULL;
	Coordinate bpCount = 0;
	const int maxline = 5000;
	char line[5000];
	IDnum sequenceCount, sequenceIndex;
	ReadSet *reads;
	short int temp_short;
	int lineLength;

	if (file != NULL)
		velvetLog("Reading read set file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	velvetLog("%li sequences found\n", (long) sequenceCount);

	reads->readCount = sequenceCount;
	
	if (reads->readCount == 0) {
		reads->sequences = NULL;
		reads->categories = NULL;
		return reads;	
	}

	reads->sequences = callocOrExit(sequenceCount, char *);
	reads->categories = callocOrExit(sequenceCount, Category);
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {

			// Reading category info
			sscanf(line, "%*[^\t]\t%*[^\t]\t%hd",
			       &temp_short);
			reads->categories[sequenceIndex + 1] = (Category) temp_short;

			if (sequenceIndex != -1)
				reads->sequences[sequenceIndex] =
				    mallocOrExit(bpCount + 1, char);
			sequenceIndex++;
			bpCount = 0;
		} if (line[0] == 'M') {;
			// Map line
		} else {
			bpCount += (Coordinate) strlen(line) - 1;

			if (sizeof(ShortLength) == sizeof(int16_t) && bpCount > SHRT_MAX) {
				velvetLog("Read %li of length %lli, longer than limit %i\n",
				       (long) sequenceIndex + 1, (long long) bpCount, SHRT_MAX);
				velvetLog("You should modify recompile with the LONGSEQUENCES option (cf. manual)\n");
				exit(1);
			}
		}
	}

	//velvetLog("Sequence %d has length %d\n", sequenceIndex, bpCount);
	reads->sequences[sequenceIndex] =
	    mallocOrExit(bpCount + 1, char);
	fclose(file);

	// Reopen file and memorize line:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file)) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				sequence[bpCount] = '\0';
			}
			sequenceIndex++;
			bpCount = 0;
			//velvetLog("Starting to read sequence %d\n",
			//       sequenceIndex);
			sequence = reads->sequences[sequenceIndex];
		} else if (line[0] == 'M') {;
			// Map line
		} else {
			lineLength = strlen(line) - 1;
			strncpy(sequence + bpCount, line, lineLength);
			bpCount += (Coordinate) lineLength;
		}
	}

	sequence[bpCount] = '\0';
	fclose(file);
	computeSecondInPair(reads);

	velvetLog("Done\n");
	return reads;

}

void logInstructions(int argc, char **argv, char *directory)
{
	int index;
	char *logFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	FILE *logFile;
	time_t date;
	char *string;

	time(&date);
	string = ctime(&date);

	strcpy(logFilename, directory);
	strcat(logFilename, "/Log");
	logFile = fopen(logFilename, "a");

	if (logFile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to %s", logFilename);

	velvetFprintf(logFile, "%s", string);

	for (index = 0; index < argc; index++)
		velvetFprintf(logFile, " %s", argv[index]);

	velvetFprintf(logFile, "\n");

	velvetFprintf(logFile, "Version %i.%i.%2.2i\n", VERSION_NUMBER,
	       RELEASE_NUMBER, UPDATE_NUMBER);
	velvetFprintf(logFile, "Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)\n");
	velvetFprintf(logFile, "This is free software; see the source for copying conditions.  There is NO\n");
	velvetFprintf(logFile, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
	velvetFprintf(logFile, "Compilation settings:\n");
	velvetFprintf(logFile, "CATEGORIES = %i\n", CATEGORIES);
	velvetFprintf(logFile, "MAXKMERLENGTH = %i\n", MAXKMERLENGTH);
#ifdef _OPENMP
	velvetFprintf(logFile, "OPENMP\n");
#endif
#ifdef LONGSEQUENCES
	velvetFprintf(logFile, "LONGSEQUENCES\n");
#endif
#ifdef BIGASSEMBLY
	velvetFprintf(logFile, "BIGASSEMBLY\n");
#endif
#ifdef COLOR
	velvetFprintf(logFile, "COLOR\n");
#endif
#ifdef DEBUG
	velvetFprintf(logFile, "DEBUG\n");
#endif
	velvetFprintf(logFile, "\n");

	fclose(logFile);
	free(logFilename);
}

void destroyReadSet(ReadSet * reads)
{
	IDnum index;

	if (reads == NULL)
		return;

	if (reads->sequences != NULL)
	{
		for (index = 0; index < reads->readCount; index++)
			free(reads->sequences[index]);
		free(reads->sequences);
	}

	if (reads->tSequences != NULL)
		free (reads->tSequences);

	if (reads->tSeqMem != NULL)
		free (reads->tSeqMem);

	if (reads->labels != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->labels[index]);

	if (reads->confidenceScores != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->confidenceScores[index]);

	if (reads->kmerProbabilities != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->kmerProbabilities[index]);

	free(reads->labels);
	free(reads->confidenceScores);
	free(reads->kmerProbabilities);
	free(reads->mateReads);
	free(reads->categories);
	free(reads->secondInPair);
	free(reads);
}

ShortLength *getSequenceLengths(ReadSet * reads, int wordLength)
{
	ShortLength *lengths = callocOrExit(reads->readCount, ShortLength);
	IDnum index;
	int lengthOffset = wordLength - 1;

	for (index = 0; index < reads->readCount; index++)
		lengths[index] =
		    getLength(getTightStringInArray(reads->tSequences, index)) - lengthOffset;

	return lengths;
}

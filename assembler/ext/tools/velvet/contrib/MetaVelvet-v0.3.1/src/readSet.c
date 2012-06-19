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

#include "globals.h"
#include "tightString.h"
#include "readSet.h"
#include "utility.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include "../../third-party/zlib-1.2.3/Win32/include/zlib.h"
#else
#include "../../third-party/zlib-1.2.3/zlib.h"
#endif

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

ReadSet *newReadSet()
{
	ReadSet *rs = callocOrExit(1, ReadSet);
	return rs;
}

static void velvetifySequence(char * str) {
	int i = strlen(str) - 1;
	char c;

	for (i = strlen(str) - 1; i >= 0; i--) {
		c = str[i];
		switch (c) {
		case '\n':
		case '\r':
		case EOF:
			str[i] = '\0';
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
			str[i] = 'A';
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
}

static void writeFastaSequence(FILE * outfile, const char * str)
{
       size_t length = strlen(str);
       size_t start;
       for (start = 0; start < length; start += 60)
	       fprintf(outfile, "%.60s\n", &str[start]);
}
 
ReadSet *newReadSetAroundTightStringArray(TightString ** array,
					  IDnum length)
{
	ReadSet *rs = newReadSet();
	rs->tSequences = array;
	rs->readCount = length;
	return rs;
}

void concatenateReadSets(ReadSet * A, ReadSet * B)
{
	ReadSet tmp;
	IDnum index;

	// Read count:  
	tmp.readCount = A->readCount + B->readCount;

	// Sequences
	if (A->sequences != NULL || B->sequences != NULL) {
		tmp.sequences = mallocOrExit(tmp.readCount, char *);
		if (A->sequences != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.sequences[index] = A->sequences[index];
			free(A->sequences);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.sequences[index] = NULL;

		if (B->sequences != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.sequences[A->readCount + index] =
				    B->sequences[index];
			free(B->sequences);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.sequences[A->readCount + index] = NULL;
	} else
		tmp.sequences = NULL;

	// tSequences
	if (A->tSequences != NULL || B->tSequences != NULL) {
		tmp.tSequences =
		    mallocOrExit(tmp.readCount, TightString *);

		if (A->tSequences != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.tSequences[index] =
				    A->tSequences[index];
			free(A->tSequences);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.tSequences[index] = NULL;

		if (B->tSequences != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.tSequences[A->readCount + index] =
				    B->tSequences[index];
			free(B->tSequences);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.tSequences[A->readCount + index] =
				    NULL;
	} else
		tmp.tSequences = NULL;

	// Labels
	if (A->labels != NULL || B->labels != NULL) {
		tmp.labels = mallocOrExit(tmp.readCount, char *);

		if (A->labels != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.labels[index] = A->labels[index];
			free(A->labels);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.labels[index] = NULL;

		if (B->labels != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.labels[A->readCount + index] =
				    B->labels[index];
			free(B->labels);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.labels[A->readCount + index] = NULL;
	} else
		tmp.labels = NULL;


	// Confidence scores
	if (A->confidenceScores != NULL || B->confidenceScores != NULL) {
		tmp.confidenceScores =
		    mallocOrExit(tmp.readCount, Quality *);

		if (A->confidenceScores != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.confidenceScores[index] =
				    A->confidenceScores[index];
			free(A->confidenceScores);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.confidenceScores[index] = NULL;

		if (B->confidenceScores != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.confidenceScores[A->readCount +
						     index] =
				    B->confidenceScores[index];
			free(B->confidenceScores);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.confidenceScores[A->readCount +
						     index] = NULL;
	} else
		tmp.confidenceScores = NULL;

	// Kmer probabilities 
	if (A->kmerProbabilities != NULL || B->kmerProbabilities != NULL) {
		tmp.kmerProbabilities =
		    mallocOrExit(tmp.readCount, Quality *);

		if (A->kmerProbabilities != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.kmerProbabilities[index] =
				    A->kmerProbabilities[index];
			free(A->kmerProbabilities);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.kmerProbabilities[index] = NULL;

		if (B->kmerProbabilities != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.kmerProbabilities[A->readCount +
						      index] =
				    B->kmerProbabilities[index];
			free(B->kmerProbabilities);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.kmerProbabilities[A->readCount +
						      index] = NULL;
	} else
		tmp.kmerProbabilities = NULL;

	// Mate reads 
	if (A->mateReads != NULL || B->mateReads != NULL) {
		tmp.mateReads = mallocOrExit(tmp.readCount, IDnum);

		if (A->mateReads != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.mateReads[index] = A->mateReads[index];
			free(A->mateReads);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.mateReads[index] = 0;

		if (B->mateReads != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.mateReads[A->readCount + index] =
				    B->mateReads[index];
			free(B->mateReads);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.mateReads[A->readCount + index] = 0;
	} else
		tmp.mateReads = NULL;

	// Categories
	if (A->categories != NULL || B->categories != NULL) {
		tmp.categories = mallocOrExit(tmp.readCount, Quality *);

		if (A->categories != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.categories[index] =
				    A->categories[index];
			free(A->categories);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.categories[index] = CATEGORIES;

		if (B->categories != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.categories[A->readCount + index] =
				    B->categories[index];
			free(B->categories);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.categories[A->readCount + index] =
				    CATEGORIES;
	} else
		tmp.categories = NULL;

	// Put everything back into A
	A->readCount = tmp.readCount;
	A->sequences = tmp.sequences;
	A->tSequences = tmp.tSequences;
	A->labels = tmp.labels;
	A->confidenceScores = tmp.confidenceScores;
	A->kmerProbabilities = tmp.kmerProbabilities;
	A->mateReads = tmp.mateReads;
	A->categories = tmp.categories;

	// Deallocate
	free(B);
}

void convertSequences(ReadSet * rs)
{
	rs->tSequences =
	    newTightStringArrayFromStringArray(rs->sequences,
					       rs->readCount);
	rs->sequences = NULL;
}

static Probability convertQualityScore(Quality score)
{
	return (Probability) 1 - pow(10, -score / ((double) 10));
}

void convertConfidenceScores(ReadSet * rs, int WORDLENGTH)
{
	Quality *baseCallerScores;
	Probability *kmerProbabilities;
	IDnum index;
	Coordinate position;
	Probability proba;

	rs->kmerProbabilities =
	    mallocOrExit(rs->readCount, Probability *);

	for (index = 0; index < rs->readCount; index++) {
		rs->kmerProbabilities[index] =
		    mallocOrExit(getLength(rs->tSequences[index]) - WORDLENGTH +
			    1, Probability);
		kmerProbabilities = rs->kmerProbabilities[index];
		baseCallerScores = rs->confidenceScores[index];

		proba = 1;
		for (position = 0;
		     position < getLength(rs->tSequences[index]);
		     position++) {
			proba *=
			    convertQualityScore(baseCallerScores
						[position]);
			if (position < WORDLENGTH)
				continue;

			proba /=
			    convertQualityScore(baseCallerScores
						[position - WORDLENGTH]);
			kmerProbabilities[position - WORDLENGTH + 1] =
			    proba;
		}

		rs->confidenceScores[index] = NULL;
		free(baseCallerScores);
	}

	free(rs->confidenceScores);
	rs->confidenceScores = NULL;
}

void categorizeReads(ReadSet * readSet, Category category)
{
	IDnum index;

	if (readSet->categories == NULL) 
		readSet->categories =
		    mallocOrExit(readSet->readCount, Category);

	for (index = 0; index < readSet->readCount; index++)
		readSet->categories[index] = category;
}

void simplifyReads(ReadSet * readSet)
{
	IDnum index;

	if (readSet->categories == NULL)
		readSet->categories =
		    mallocOrExit(readSet->readCount, Category);

	for (index = 0; index < readSet->readCount; index++) {
		if (readSet->categories[index] < CATEGORIES) {
			readSet->categories[index] /= 2;
			readSet->categories[index] *= 2;
		}
	}
}

void exportIDMapping(char *filename, ReadSet * reads)
{
	IDnum index;
	FILE *outfile = fopen(filename, "w");

	if (outfile == NULL) {
		printf("Couldn't open %s, sorry\n", filename);
		return;
	} else
		puts("Writing into file...");

	if (reads->labels == NULL) {
		fclose(outfile);
		return;
	}

	for (index = 0; index < reads->readCount; index++)
		if (reads->labels != NULL)
			fprintf(outfile, "s/SEQUENCE %ld/%s/\n", (long) (index + 1),
				reads->labels[index]);

	fclose(outfile);

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

// Imports sequences from a fastq file 
// Memory space allocated within this function.
static void readSolexaFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	FILE *file = fopen(filename, "r");
	IDnum counter = 0;
	const int maxline = 500;
	char line[500];
	char readName[500];
	char readSeq[500];
	char str[100];
	Coordinate start;

	if (strcmp(filename, "-"))
		file = fopen(filename, "r");
	else
		file = stdin;

	if (file != NULL)
		printf("Reading Solexa file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	while (fgets(line, maxline, file) != NULL)
		if (strchr(line, '.') == NULL) {
			sscanf(line, "%s\t%*i\t%*i\t%*i\t%*c%[^\n]",
			       readName, readSeq);
			fprintf(outfile, ">%s\t%ld\t%d\n", readName, (long) ((*sequenceIndex)++), (int) cat);
			velvetifySequence(readSeq);
			start = 0;
			while (start <= strlen(readSeq)) {
				strncpy(str, readSeq + start, 60);
				str[60] = '\0';
				fprintf(outfile, "%s\n", str);
				start += 60;
			}

			counter++;
		}

	fclose(file);

	printf("%d sequences found\n", counter);
	puts("Done");
}

static void readElandFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	FILE *file = fopen(filename, "r");
	IDnum counter = 0;
	const int maxline = 5000;
	char line[5000];
	char readName[5000];
	char readSeq[5000];
	char str[100];
	Coordinate start;

	if (strcmp(filename, "-"))
		file = fopen(filename, "r");
	else
		file = stdin;

	if (file != NULL)
		printf("Reading Solexa file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Reopen file and memorize line:
	while (fgets(line, maxline, file) != NULL) {
		sscanf(line, "%[^\t]\t%[^\t\n]",
		       readName, readSeq);
		fprintf(outfile, ">%s\t%ld\t%d\n", readName, (long) ((*sequenceIndex)++), (int) cat);
		velvetifySequence(readSeq);
		start = 0;
		while (start <= strlen(readSeq)) {
			strncpy(str, readSeq + start, 60);
			str[60] = '\0';
			fprintf(outfile, "%s\n", str);
			start += 60;
		}

		counter++;
	}

	fclose(file);

	printf("%d sequences found\n", counter);
	puts("Done");
}

void goToEndOfLine(char *line, FILE * file)
{
	size_t length = strlen(line);
	char c = line[length - 1];

	while (c != '\n')
		c = fgetc(file);
}

// Imports sequences from a fastq file 
// Memory space allocated within this function.
static void readFastQFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	FILE *file;
	const int maxline = 5000;
	char line[5000];
	char str[100];
	IDnum counter = 0;
	Coordinate start, i;
	char c;

	if (strcmp(filename, "-"))
		file = fopen(filename, "r");
	else 
		file = stdin;

	if (file != NULL)
		printf("Reading FastQ file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Checking if FastQ
	c = getc(file);
	if (c != '@') 
		exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastQ format", filename);
	ungetc(c, file);	

	while(fgets(line, maxline, file)) { 

		for (i = strlen(line) - 1;
		     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
			line[i] = '\0';
		}

		fprintf(outfile,">%s\t%ld\t%d\n", line + 1, (long) ((*sequenceIndex)++), (int) cat);
		counter++;

		if(!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", filename);

		velvetifySequence(line);

		start = 0;
		while (start <= strlen(line)) {
			strncpy(str, line + start, 60);
			str[60] = '\0';
			fprintf(outfile, "%s\n", str);
			start += 60;
		}

		if(!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", filename);
		if(!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", filename);
	}

	fclose(file);
	printf("%d reads found.\n", counter);
	puts("Done");
}

// Imports sequences from a zipped rfastq file 
// Memory space allocated within this function.
static void readFastQGZFile(FILE * outfile, char *filename, Category cat, IDnum *sequenceIndex)
{
	gzFile file;
	const int maxline = 5000;
	char line[5000];
	char str[100];
	IDnum counter = 0;
	Coordinate start, i;
	char c;

	if (strcmp(filename, "-"))
		file = gzopen(filename, "rb");
	else { 
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		printf("Reading FastQ file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Checking if FastQ
	c = gzgetc(file);
	if (c != '@') 
		exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastQ format", filename);
	gzungetc(c, file);	

	while (gzgets(file, line, maxline)) {
		for (i = strlen(line) - 1;
		     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
			line[i] = '\0';
		}

		fprintf(outfile,">%s\t%ld\t%d\n", line + 1, (long) ((*sequenceIndex)++), (int) cat);
		counter++;

		gzgets(file, line, maxline);

		velvetifySequence(line);

		start = 0;
		while (start <= strlen(line)) {
			strncpy(str, line + start, 60);
			str[60] = '\0';
			fprintf(outfile, "%s\n", str);
			start += 60;
		}

		gzgets(file, line, maxline);
		gzgets(file, line, maxline);
	}

	gzclose(file);
	printf("%d reads found.\n", counter);
	puts("Done");
}

// Imports sequences from a fasta file 
// Memory is allocated within the function 
static void readFastAFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	FILE *file;
	const int maxline = 5000;
	char line[5000];
	char str[100];
	IDnum counter = 0;
	Coordinate i;
	char c;
	Coordinate start;
	int offset = 0;

	if (strcmp(filename, "-"))
		file = fopen(filename, "r");
	else
		file = stdin;

	if (file != NULL)
		printf("Reading FastA file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Checking if FastA
	c = getc(file);
	if (c != '>') 
		exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastA format", filename);
	ungetc(c, file);	

	while (fgets(line, maxline, file)) {
		if (line[0] == '>') {
			if (strchr(line,'\t')) {
				printf("FastA headers in %s contain tabs, please remove them.\n", filename);
				printf("E.g.: %s", line);
				puts("Exiting");
				exit(1);
			}

			if (offset != 0) { 
				fprintf(outfile, "\n");
				offset = 0;
			}

			for (i = strlen(line) - 1;
			     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
				line[i] = '\0';
			}

			fprintf(outfile,"%s\t%ld\t%d\n", line, (long) ((*sequenceIndex)++), (int) cat);
			counter++;
		} else {
			velvetifySequence(line);
			start = 0;
			while (start < strlen(line)) {
				strncpy(str, line + start, 60 - offset);
				str[60 - offset] = '\0';
				fprintf(outfile, "%s", str);
				offset += strlen(str);
				if (offset >= 60) {
					fprintf(outfile, "\n");
					offset = 0;
				}
				start += strlen(str);
			}
		}
	}

	if (offset != 0) 
		fprintf(outfile, "\n");
	fclose(file);

	printf("%d sequences found\n", counter);
	puts("Done");
}

// Imports sequences from a zipped fasta file 
// Memory is allocated within the function 
static void readFastAGZFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	gzFile file;
	const int maxline = 5000;
	char line[5000];
	char str[100];
	IDnum counter = 0;
	Coordinate i, start;
	char c;
	int offset = 0;

	if (strcmp(filename, "-"))
		file = gzopen(filename, "rb");
	else { 
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		printf("Reading zipped FastA file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Checking if FastA
	c = gzgetc(file);
	if (c != '>') 
		exitErrorf(EXIT_FAILURE, false, "%s does not seem to be in FastA format", filename);
	gzungetc(c, file);	

	while (gzgets(file, line, maxline)) {
		if (line[0] == '>') {
			if (offset != 0) { 
				fprintf(outfile, "\n");
				offset = 0;
			}

			for (i = strlen(line) - 1;
			     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
				line[i] = '\0';
			}

			fprintf(outfile, "%s\t%ld\t%d\n", line, (long) ((*sequenceIndex)++), (int) cat);	
			counter++;
		} else {
			velvetifySequence(line);

			start = 0;
			while (start < strlen(line)) {
				strncpy(str, line + start, 60 - offset);
				str[60 - offset] = '\0';
				fprintf(outfile, "%s", str);
				offset += strlen(str);
				if (offset >= 60) {
					fprintf(outfile, "\n");
					offset = 0;
				}
				start += strlen(str);
			}
		}
	}

	if (offset != 0) 
		fprintf(outfile, "\n");
	gzclose(file);

	printf("%d sequences found\n", counter);
	puts("Done");
}

// Parser for new output
static void readMAQGZFile(FILE* outfile, char *filename, Category cat, IDnum * sequenceIndex)
{
	gzFile file;
	const int maxline = 1000;
	char line[1000];
	IDnum counter = 0;
	char readName[500];
	char readSeq[500];
	char str[100];
	Coordinate start;

	if (strcmp(filename, "-"))
		file = gzopen(filename, "rb");
	else { 
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		printf("Reading zipped MAQ file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Reopen file and memorize line:
	while (gzgets(file, line, maxline)) {
		sscanf(line, "%s\t%*i\t%*i\t%*c\t%*i\t%*i\t%*i\t%*i\t%*i\t%*i\t%*i\t%*i\t%*i\t%*i\t%[^\t]",
		       readName, readSeq);
		fprintf(outfile, ">%s\t%ld\t%d\n", readName, (long) ((*sequenceIndex)++), (int) cat);
		velvetifySequence(readSeq);
		start = 0;
		while (start <= strlen(readSeq)) {
			strncpy(str, readSeq + start, 60);
			str[60] = '\0';
			fprintf(outfile, "%s\n", str);
			start += 60;
		}

		counter++;
	}

	gzclose(file);

	printf("%d sequences found\n", counter);
	puts("Done");
}

static void readSAMFile(FILE *outfile, char *filename, Category cat, IDnum *sequenceIndex)
{
	char line[5000];
	unsigned long lineno, readCount;
	char previous_qname_pairing[10];
	char previous_qname[5000];
	char previous_seq[5000];
	boolean previous_paired = false;

	FILE *file = (strcmp(filename, "-") != 0)? fopen(filename, "r") : stdin;
	if (file)
		printf("Reading SAM file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	readCount = 0;
	for (lineno = 1; fgets(line, sizeof(line), file); lineno++)
		if (line[0] != '@') {
			char *qname, *flag, *seq;
			int i;

			qname = strtok(line, "\t");
			flag  = strtok(NULL, "\t");
			for (i = 3; i < 10; i++)
				(void) strtok(NULL, "\t");
			seq = strtok(NULL, "\t");

			if (seq == NULL) {
				fprintf(stderr,
					"Line #%lu: ignoring SAM record with too few fields\n",
					lineno);
			}
			else if (strcmp(seq, "*") == 0) {
				fprintf(stderr,
					"Line #%lu: ignoring SAM record with omitted SEQ field\n",
					lineno);
			}
			else {
				// Accept flags represented in either decimal or hex:
				int flagbits = strtol(flag, NULL, 0);

				const char *qname_pairing = "";
				if (flagbits & 0x40)
					qname_pairing = "/1";
				else if (flagbits & 0x80)
					qname_pairing = "/2";

				if (flagbits & 0x10)
					reverseComplementSequence(seq);

				// Determine if paired to previous read
				if (readCount > 0) {
					if (cat % 2) {
						if (previous_paired) {
							// Last read paired to penultimate read
							fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
								(long) ((*sequenceIndex)++), (int) cat);
							writeFastaSequence(outfile, previous_seq);
							previous_paired = false;
						} else if (strcmp(qname, previous_qname) == 0) {
							// Last read paired to current reads
							fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
								(long) ((*sequenceIndex)++), (int) cat);
							writeFastaSequence(outfile, previous_seq);
							previous_paired = true;
						} else {
							// Last read unpaired
							fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
								(long) ((*sequenceIndex)++), (int) cat - 1);
							writeFastaSequence(outfile, previous_seq);
							previous_paired = false;
						}
					} else {
						// Unpaired dataset 
						fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
							(long) ((*sequenceIndex)++), (int) cat);
						writeFastaSequence(outfile, previous_seq);
					}
				}

				strcpy(previous_qname, qname);
				strcpy(previous_qname_pairing, qname_pairing);
				strcpy(previous_seq, seq);
				velvetifySequence(previous_seq);

				readCount++;
			}
		}

	if (readCount) {
		if (cat % 2) {
			if (previous_paired) {
				// Last read paired to penultimate read
				fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
					(long) ((*sequenceIndex)++), (int) cat);
				writeFastaSequence(outfile, previous_seq);
			} else {
				// Last read unpaired
				fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
					(long) ((*sequenceIndex)++), (int) cat - 1);
				writeFastaSequence(outfile, previous_seq);
			}
		} else {
			// Unpaired dataset 
			fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
				(long) ((*sequenceIndex)++), (int) cat);
			writeFastaSequence(outfile, previous_seq);
		}

	}

	fclose(file);
	printf("%lu reads found.\nDone\n", readCount);
}

static int readBAMint32(gzFile file)
{
	unsigned char buffer[4];
	if (gzread(file, buffer, 4) != 4)
		exitErrorf(EXIT_FAILURE, false, "BAM file header truncated");

	return int32(buffer);
}

static void readBAMFile(FILE *outfile, char *filename, Category cat, IDnum *sequenceIndex)
{
	size_t seqCapacity = 0;
	char *seq = NULL;
	size_t bufferCapacity = 4;
	unsigned char *buffer = mallocOrExit(bufferCapacity, unsigned char);
	unsigned long recno, readCount;
	int i, refCount;
	gzFile file;
	char previous_qname_pairing[10];
	char previous_qname[5000];
	char previous_seq[5000];
	boolean previous_paired = false;

	if (strcmp(filename, "-") != 0)
		file = gzopen(filename, "rb");
	else {
		file = gzdopen(fileno(stdin), "rb");
		SET_BINARY_MODE(stdin);
	}

	if (file != NULL)
		printf("Reading BAM file %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	if (! (gzread(file, buffer, 4) == 4 && memcmp(buffer, "BAM\1", 4) == 0))
		exitErrorf(EXIT_FAILURE, false, "%s is not in BAM format", filename);

	// Skip header text
	if (gzseek(file, readBAMint32(file), SEEK_CUR) == -1)
		exitErrorf(EXIT_FAILURE, false, "gzseek failed");

	// Skip header reference list
	refCount = readBAMint32(file);
	for (i = 0; i < refCount; i++) {
		if (gzseek(file, readBAMint32(file) + 4, SEEK_CUR) == -1)
			exitErrorf(EXIT_FAILURE, false, "gzseek failed");
	}

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
			fprintf(stderr,
				"Record #%lu: ignoring BAM record with omitted SEQ field\n",
				recno);
		}
		else {
			int readNameLength = buffer[12];
			int flag_nc = int32(&buffer[16]);
			int flagbits = flag_nc >> 16;
			int cigarLength = flag_nc & 0xffff;
			char *qname = (char *)&buffer[36];
			unsigned char *rawseq =
					&buffer[36 + readNameLength + 4 * cigarLength];

			const char *qname_pairing = "";
			if (flagbits & 0x40)
				qname_pairing = "/1";
			else if (flagbits & 0x80)
				qname_pairing = "/2";

			if (seqCapacity < readLength + 1) {
				seqCapacity = readLength * 2 + 1;
				seq = reallocOrExit(seq, seqCapacity, char);
			}

			for (i = 0; i < readLength; i += 2) {
				static const char decode_bases[] = "=ACMGRSVTWYHKDBN";
				unsigned int packed = *rawseq++;
				seq[i] = decode_bases[packed >> 4];
				seq[i+1] = decode_bases[packed & 0xf];
			}
			seq[readLength] = '\0';

			if (flagbits & 0x10)
				reverseComplementSequence(seq);

			// Determine if paired to previous read
			if (readCount > 0) {
				if (cat % 2) {
					 if (previous_paired) {
						// Last read paired to penultimate read
						fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
							(long) ((*sequenceIndex)++), (int) cat);
						writeFastaSequence(outfile, previous_seq);
						previous_paired = false;
					} else if (strcmp(qname, previous_qname) == 0) {
						// Last read paired to current reads
						fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
							(long) ((*sequenceIndex)++), (int) cat);
						writeFastaSequence(outfile, previous_seq);
						previous_paired = true;
					} else {
						// Last read unpaired
						fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
							(long) ((*sequenceIndex)++), (int) cat - 1);
						writeFastaSequence(outfile, previous_seq);
						previous_paired = false;
					}
				} else {
					// Unpaired dataset 
					fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
						(long) ((*sequenceIndex)++), (int) cat);
					writeFastaSequence(outfile, previous_seq);
				}
			}

			strcpy(previous_qname, qname);
			strcpy(previous_qname_pairing, qname_pairing);
			strcpy(previous_seq, seq);
			velvetifySequence(previous_seq);

			readCount++;
		}
	}

	if (readCount) {
		if (cat % 2) {
			if (previous_paired) {
				// Last read paired to penultimate read
				fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
					(long) ((*sequenceIndex)++), (int) cat);
				writeFastaSequence(outfile, previous_seq);
			} else {
				// Last read unpaired
				fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
					(long) ((*sequenceIndex)++), (int) cat - 1);
				writeFastaSequence(outfile, previous_seq);
			}
		} else {
			// Unpaired dataset 
			fprintf(outfile, ">%s%s\t%ld\t%d\n", previous_qname, previous_qname_pairing,
				(long) ((*sequenceIndex)++), (int) cat);
			writeFastaSequence(outfile, previous_seq);
		}
	}

	free(seq);
	free(buffer);

	gzclose(file);
	printf("%lu reads found.\nDone\n", readCount);
}

#define FASTQ 1
#define FASTA 2
#define GERALD 3
#define ELAND 4
#define FASTA_GZ 5
#define FASTQ_GZ 6
#define MAQ_GZ 7
#define SAM 8
#define BAM 9

// General argument parser for most functions
// Basically a reused portion of toplevel code dumped into here
void parseDataAndReadFiles(char * filename, int argc, char **argv, boolean * double_strand)
{
	int argIndex = 1;
	FILE *outfile = fopen(filename, "w");
	int filetype = FASTA;
	Category cat = 0;
	IDnum sequenceIndex = 1;
	short short_var;

	if (argc < 2) {
		puts("Wrong number of arguments!");
		puts("Correct usage:");
		puts("run -<filetype> <list of files> [...] ");
		puts("Allowed filetypes:");
		puts("\t-fasta");
		puts("\t-fastq");
		puts("\t-solexa");
		puts("\t-eland");
		puts("If reading exclusively fasta file, the -fasta parameter is not necessary");
		exit(1);
	}

	for (argIndex = 1; argIndex < argc; argIndex++) {
		if (argv[argIndex][0] == '-' && strlen(argv[argIndex]) > 1) {

			if (strcmp(argv[argIndex], "-fastq") == 0)
				filetype = FASTQ;
			else if (strcmp(argv[argIndex], "-fasta") == 0)
				filetype = FASTA;
			else if (strcmp(argv[argIndex], "-gerald") == 0)
				filetype = GERALD;
			else if (strcmp(argv[argIndex], "-eland") == 0)
				filetype = ELAND;
			else if (strcmp(argv[argIndex], "-fastq.gz") == 0)
				filetype = FASTQ_GZ;
			else if (strcmp(argv[argIndex], "-fasta.gz") == 0)
				filetype = FASTA_GZ;
			else if (strcmp(argv[argIndex], "-sam") == 0)
				filetype = SAM;
			else if (strcmp(argv[argIndex], "-bam") == 0)
				filetype = BAM;
			else if (strcmp(argv[argIndex], "-maq.gz") == 0)
				filetype = MAQ_GZ;
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
					printf("Unknown option: %s\n",
					       argv[argIndex]);
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
					printf("Unknown option: %s\n",
					       argv[argIndex]);
					exit(1);
				}
				cat--;
				cat *= 2;
			} else if (strcmp(argv[argIndex], "-long") == 0)
				cat = CATEGORIES * 2;
			else if (strcmp(argv[argIndex], "-longPaired") ==
				 0)
				cat = CATEGORIES * 2 + 1;
			else if (strcmp(argv[argIndex], "-strand_specific") 
				 == 0)
				*double_strand = false;
			else {
				printf("Unknown option: %s\n",
				       argv[argIndex]);
				exit(1);
			}

			continue;
		}

		if (cat == -1)
			continue;

		switch (filetype) {
		case FASTA:
			readFastAFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case FASTQ:
			readFastQFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case GERALD:
			readSolexaFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case ELAND:
			readElandFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case FASTA_GZ:
			readFastAGZFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case FASTQ_GZ:
			readFastQGZFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case SAM:
			readSAMFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case BAM:
			readBAMFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		case MAQ_GZ:
			readMAQGZFile(outfile, argv[argIndex], cat, &sequenceIndex);
			break;
		default:
			puts("Screw up in parser... exiting");
			exit(1);
		}
	}

	fclose(outfile);
}

void createReadPairingArray(ReadSet* reads) {
	IDnum index;
	IDnum *mateReads = mallocOrExit(reads->readCount, IDnum);

	for (index = 0; index < reads->readCount; index++) 
		mateReads[index] = -1;

	reads->mateReads = mateReads;
}

boolean pairUpReads(ReadSet * reads, Category cat)
{
	int phase = 0;
	IDnum index;
	boolean found = false;

	for (index = 0; index < reads->readCount; index++) {
		if (reads->categories[index] != cat) {
			if (phase == 1) {
				reads->mateReads[index - 1] = -1;
				reads->categories[index - 1]--;
				phase = 0;
			}
		} else if (phase == 0) {
			found = true;
			reads->mateReads[index] = index + 1;
			phase = 1;
		} else {
			found = true;
			reads->mateReads[index] = index - 1;
			phase = 0;
		}
	}

	return found;
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
		if (!dubiousReads[index])
			continue;

		pairID = mateReads[index];

		if (pairID != -1) {
			//printf("Separating %d and %d\n", index, pairID);
			mateReads[index] = -1;
			mateReads[pairID] = -1;
		}
	}
}

static void exportRead(FILE * outfile, ReadSet * reads, IDnum index)
{
	Coordinate start, finish;
	char str[100];
	TightString *sequence = reads->tSequences[index];

	if (sequence == NULL)
		return;

	fprintf(outfile, ">SEQUENCE_%ld_length_%lld", (long) index,
		(long long) getLength(sequence));

	if (reads->categories != NULL)
		fprintf(outfile, "\t%i", (int) reads->categories[index]);

	fprintf(outfile, "\n");

	start = 0;
	while (start <= getLength(sequence)) {
		finish = start + 60;
		readTightStringFragment(sequence, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}

	fflush(outfile);
}

void exportReadSet(char *filename, ReadSet * reads)
{
	IDnum index;
	FILE *outfile = fopen(filename, "w+");

	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		printf("Writing into readset file: %s\n", filename);

	for (index = 0; index < reads->readCount; index++) {
		exportRead(outfile, reads, index);
	}

	fclose(outfile);

	puts("Done");
}

ReadSet *importReadSet(char *filename)
{
	FILE *file = fopen(filename, "r");
	char *sequence = NULL;
	Coordinate bpCount = 0;
	const int maxline = 5000;
	char line[5000];
	IDnum sequenceCount, sequenceIndex;
	IDnum index;
	ReadSet *reads;
	short int temp_short;

	if (file != NULL)
		printf("Reading read set file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	printf("%d sequences found\n", sequenceCount);

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
		} else {
			bpCount += (Coordinate) strlen(line) - 1;
		}
	}

	//printf("Sequence %d has length %d\n", sequenceIndex, bpCount);
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
			//printf("Starting to read sequence %d\n",
			//       sequenceIndex);
			sequence = reads->sequences[sequenceIndex];
		} else {
			for (index = 0; index < (Coordinate) strlen(line) - 1;
			     index++)
				sequence[bpCount + index] = line[index];
			bpCount += (Coordinate) (strlen(line) - 1);
		}
	}

	sequence[bpCount] = '\0';
	fclose(file);

	puts("Done");
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

	fprintf(logFile, "%s", string);

	for (index = 0; index < argc; index++)
		fprintf(logFile, " %s", argv[index]);

	fprintf(logFile, "\n");

	fclose(logFile);
	free(logFilename);
}

void destroyReadSet(ReadSet * reads)
{
	IDnum index;

	if (reads == NULL)
		return;

	if (reads->sequences != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->sequences[index]);

	if (reads->tSequences != NULL)
		for (index = 0; index < reads->readCount; index++)
			destroyTightString(reads->tSequences[index]);

	if (reads->labels != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->labels[index]);

	if (reads->confidenceScores != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->confidenceScores[index]);

	if (reads->kmerProbabilities != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->kmerProbabilities[index]);

	free(reads->sequences);
	free(reads->tSequences);
	free(reads->labels);
	free(reads->confidenceScores);
	free(reads->kmerProbabilities);
	free(reads->mateReads);
	free(reads->categories);
	free(reads);
}

Coordinate *getSequenceLengths(ReadSet * reads, int wordLength)
{
	Coordinate *lengths = callocOrExit(reads->readCount, Coordinate);
	IDnum index;
	int lengthOffset = wordLength - 1;

	for (index = 0; index < reads->readCount; index++)
		lengths[index] =
		    getLength(reads->tSequences[index]) - lengthOffset;

	return lengths;
}

Coordinate *getSequenceLengthsFromFile(char *filename, int wordLength)
{
	Coordinate *lengths;
	FILE *file = fopen(filename, "r");
	Coordinate bpCount = 0;
	const int maxline = 100;
	char line[100];
	IDnum sequenceCount, sequenceIndex;
	int lengthOffset = wordLength - 1;

	if (file != NULL)
		printf("Reading read set file %s;\n", filename);
	else 
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);

	lengths = callocOrExit(sequenceCount, Coordinate);
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {
			if (sequenceIndex != -1)
				lengths[sequenceIndex] =
				    bpCount - lengthOffset;
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (Coordinate) strlen(line) - 1;
		}
	}
	lengths[sequenceIndex] =  bpCount - lengthOffset;
	fclose(file);

	return lengths;
}

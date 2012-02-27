/***************************************************************************
 * Title:          SFF2Fasta.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>
#include "utils.h"

using namespace std;

 
uint16_t BigEndian16(uint16_t n)
{
	unsigned char* buf = (unsigned char*) &n;
	return (uint16_t)
		(((uint16_t)buf[1]) +
		 ((uint16_t)buf[0]<<8));
}
  
uint32_t BigEndian32(uint32_t n) {
	unsigned char* buf = (unsigned char*) &n;
	
	return (uint32_t)
		(((uint32_t)buf[3]) +
		 ((uint32_t)buf[2]<<8) +
		 ((uint32_t)buf[1]<<16) +
		 ((uint32_t)buf[0]<<24));
}

uint64_t BigEndian64(uint64_t n){ 
	unsigned char* buf = (unsigned char*) &n;

	return (uint64_t) 
		(((uint64_t) buf[7] ) +
		 ((uint64_t) buf[6] << 8 ) +
		 ((uint64_t) buf[5] << 16) +
		 ((uint64_t) buf[4] << 24) +
		 ((uint64_t) buf[3] << 32) +
		 ((uint64_t) buf[2] << 40) +
		 ((uint64_t) buf[1] << 48) +
		 ((uint64_t) buf[0] << 56));
}

void SkipPadding(std::ifstream &readsIn, int numToSkip) {
	// discard the word alignment padding.
	unsigned char byte;
	int i;
	for (i = 0; i < numToSkip; i++ ){
		readsIn.read((char*) &byte, sizeof(char));
	}
}


void ReadReadHeader(std::ifstream &readsIn, 
										char* &nameP, int &nameLengthP,
										int &numBasesP, 
										int &clipQualLeftP, int &clipQualRightP,
										int &clipAdapLeftP, int &clipAdapRightP) {
	/*
		
		# read_header_length uint16_t
		# name_length uint16_t
		# number_of_bases uint32_t
		# clip_qual_left uint16_t
		# clip_qual_right uint16_t
		# clip_adapter_left uint16_t
		# clip_adapter_right uint16_t
		# name char[name_length]
		# eight_byte_padding uint8_t[*]
	*/
	uint16_t headerLength;
	uint16_t nameLength;
	uint32_t numBases;
	uint16_t clipQualLeft;
	uint16_t clipQualRight;
	uint16_t clipAdapterLeft;
	uint16_t clipAdapterRight;
	char *name;
	
	readsIn.read((char*) &headerLength, sizeof(uint16_t));
	readsIn.read((char*) &nameLength, sizeof(uint16_t));
	readsIn.read((char*) &numBases, sizeof(uint32_t));
	readsIn.read((char*) &clipQualLeft, sizeof(uint16_t));
	readsIn.read((char*) &clipQualRight, sizeof(uint16_t));
	readsIn.read((char*) &clipAdapterLeft, sizeof(uint16_t));
	readsIn.read((char*) &clipAdapterRight, sizeof(uint16_t));
	nameLength = BigEndian16(nameLength);
	name = new char[nameLength];
	readsIn.read((char*) name, sizeof(char)*nameLength);
	
	headerLength = BigEndian16(headerLength);
	numBases = BigEndian32(numBases);
	clipQualLeft = BigEndian16(clipQualLeft);
	clipQualRight = BigEndian16(clipQualRight);
	clipAdapterLeft = BigEndian16(clipAdapterLeft);
	clipAdapterRight = BigEndian16(clipAdapterRight);


	int numToSkip = (sizeof(uint16_t) 
									 + sizeof(uint16_t)
									 + sizeof(uint32_t)
									 + sizeof(uint16_t)
									 + sizeof(uint16_t)
									 + sizeof(uint16_t)
									 + sizeof(uint16_t));
	numToSkip += nameLength;

	numToSkip %= 8;
	if (numToSkip != 0)
		SkipPadding(readsIn, 8 - numToSkip);

	nameP = name;
	nameLengthP = nameLength;
	numBasesP = numBases;
	clipQualLeftP = clipQualLeft;
	clipQualRight = clipQualRight;
	clipAdapLeftP = clipAdapterLeft;
	clipAdapRightP = clipAdapterRight;
}

void ReadReadValues(std::ifstream &readsIn, 
										int numberOfFlows, int numberOfBases, 
										uint16_t *&flowValues, 
										uint8_t *&flowIndices,
										char *&bases,
										uint8_t *&qualScores) {
	flowValues = new uint16_t[numberOfFlows];
	flowIndices = new uint8_t[numberOfBases];
	bases = new char[numberOfBases];
	qualScores = new uint8_t[numberOfBases];

	readsIn.read((char*) flowValues,  sizeof(uint16_t)*numberOfFlows);
	readsIn.read((char*) flowIndices, sizeof(uint8_t)*numberOfBases);
	readsIn.read((char*) bases,       sizeof(char)*numberOfBases);
	readsIn.read((char*) qualScores,  sizeof(char)*numberOfBases);

	int numToSkip = (sizeof(uint16_t)*numberOfFlows
									 + sizeof(uint8_t)*numberOfBases
									 + sizeof(char)*numberOfBases
									 + sizeof(char)*numberOfBases) % 8;
	if (numToSkip != 0)
		SkipPadding(readsIn, 8-numToSkip);
}
										

int main(int argc, char* argv[]) {
	
	string sffFileName, fastaFileName;
	
	if (argc < 3) {
		cout << "usage: sff2fastq sffFile fastqFile" << endl;
		cout << "       This reads in a 454 sff file that contains both base and " << endl
				 << "       quality values and outputs the fastq file.  For every sequence" <<endl
				 << "       the fastq file 4 lines in the format: "<<endl << endl
				 << "@TITLE" <<endl
				 << "sequence" <<endl
				 << "+TITLE" << endl
				 << "quality_string" << endl << endl
				 << "  TITLE is the FASTA read title, and is repeated for both the sequence" << endl
				 << "           and the quality values.  " <<endl
				 << "  sequence is the nucleotide sequence of the read." << endl
				 << "  quality_string is the ascii-representation of the quality values, and shoudl"<<endl
				 << "           be the same length as the sequence."<<endl;
		exit(1);
	}

	sffFileName = argv[1];
	fastaFileName = argv[2];

	ifstream sffFile;
	ofstream fastaFile;
	openck(sffFileName, sffFile, std::ios::in | std::ios::binary);
	openck(fastaFileName, fastaFile, std::ios::out);


	/*
		Parse the sff file.  From the NCBI web page:

    * magic_number uint32_t
    * version char[4]
    * index_offset uint64_t
    * index_length uint32_t
    * number_of_reads uint32_t
    * header_length uint16_t
    * key_length uint16_t
    * number_of_flows_per_read uint16_t
    * flowgram_format_code uint8_t
    * flow_chars char[number_of_flows_per_read]
    * key_sequence char[key_length]
    * eight_byte_padding uint8_t[*]


	*/
	uint32_t magicNumber;
	uint64_t indexOffset;
	uint32_t indexLength;
	uint32_t numberOfReads;
	uint16_t headerLength;
	uint16_t keyLength;
	uint16_t flowsPerRead;
	uint8_t  formatCode;
	char *flowChars;
	char *keySequence;
	//UNUSED//	char *padding;
	uint32_t version;


	sffFile.read((char*) &magicNumber, sizeof(uint32_t));  //4,4
	sffFile.read((char*) &version, sizeof(uint32_t)); // 4, 8
	sffFile.read((char*) &indexOffset, sizeof(uint64_t));//8,16
	sffFile.read((char*) &indexLength, sizeof(uint32_t));//4,20
	sffFile.read((char*) &numberOfReads, sizeof(uint32_t));//4,24
	sffFile.read((char*) &headerLength, sizeof(uint16_t));//2,26
	sffFile.read((char*) &keyLength, sizeof(uint16_t));//2,28
	sffFile.read((char*) &flowsPerRead, sizeof(uint16_t));//2,30
	sffFile.read((char*) &formatCode, sizeof(uint8_t));//1,31

	magicNumber = BigEndian32(magicNumber);
	version     = BigEndian32(version);
	indexOffset = BigEndian64(indexOffset);
	indexLength = BigEndian32(indexLength);
	numberOfReads=BigEndian32(numberOfReads);
	headerLength= BigEndian16(headerLength);
	keyLength   = BigEndian16(keyLength);
	flowsPerRead= BigEndian16(flowsPerRead);
	

	flowChars = new char[flowsPerRead+1];
	flowChars[flowsPerRead] = '\0';
	sffFile.read((char*) flowChars, sizeof(char)*flowsPerRead);
	
	keySequence = new char[keyLength+1];
	keySequence[keyLength] = '\0';
	sffFile.read((char*) keySequence, sizeof(char) *keyLength);

	cout << "magic number: " << magicNumber << endl;
	cout << "version: " << version << endl;
	cout << "index offset: " << indexOffset << endl;
	cout << "index length: " << indexLength << endl;
	cout << "number of reads: " << numberOfReads << endl;
	cout << "header length: " << headerLength << endl;
	cout << "key length: " << keyLength << endl;
	cout << "flowsPerRead: " << flowsPerRead << endl;
	cout << "format code: " << formatCode << endl;
	cout << "flow chars: " << flowChars << endl;
	cout << "key sequence: " << keySequence << endl;

	int numToSkip = (31 + flowsPerRead + keyLength) % 8;
	if (numToSkip != 0) 
		SkipPadding(sffFile, 8- numToSkip);
	//
	// Read the sff file header.
	//
	
	int r;
	for (r = 0; r < numberOfReads; r++ ){ 
		char *name, *seq;
		int nameLength, seqLength;
		int clipQLeft, clipQRight, clipALeft, clipARight;
		ReadReadHeader(sffFile, name, nameLength,
									 seqLength,
									 clipQLeft, clipQRight,
									 clipALeft, clipARight);
		uint16_t *flowValues;
		uint8_t *flowIndices;
		uint8_t  *qualScores;

		ReadReadValues(sffFile, flowsPerRead, seqLength,
									 flowValues, flowIndices, seq, qualScores);

		
		string seqStr, nameStr;
		nameStr.assign(name, nameLength);
		seqStr.assign(seq, seqLength);
		fastaFile << "@" << name << endl;
		fastaFile << seqStr << endl;
		int i;
		for (i = 0; i < seqLength; i++ ) {
			//			cout << (int) qualScores[i] << " ";
			qualScores[i] = ((qualScores[i] <= 93) ? qualScores[i] : 93) + 33;
			//			cout << (int) qualScores[i] << ", ";
		}
		//		cout << endl;
		seqStr.assign((char*) qualScores, seqLength);
		fastaFile << "+" << name << endl;
		fastaFile << seqStr << endl;
		
		delete [] flowValues;
		delete [] flowIndices;
		delete [] qualScores;
		delete [] seq;
	}
	return 0;
}

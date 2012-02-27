/***************************************************************************
 * Title:          ClipReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/17/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
void ReadQualityScores(std::ifstream &scoresIn, 
											 std::string &title,
											 std::vector<ssize_t> &scores);

void PrintQualityScores(std::ofstream &scoresOut,
												std::string &title,
												std::vector<ssize_t> &scores);

int main(int argc, char* argv[]) {
	std::string seqFileName, seqOutFileName, clipFileName;
	std::string scoresFileName, scoresOutFileName;
	int argi = 1;
	if (argc != 6) {
		std::cout << "usage: clipReads seqIn qualIn clipFile seqOut qualOut " << std::endl;
		exit(1);
	}
	seqFileName       = argv[argi++];
	scoresFileName    = argv[argi++];
	clipFileName      = argv[argi++];
	seqOutFileName    = argv[argi++];
	scoresOutFileName = argv[argi++];

	std::ifstream seqFile, clipFile, scoresFile;
	std::ofstream seqOutFile, scoresOutFile;
	
	openck(seqFileName,   seqFile, std::ios::in);
	openck(scoresFileName, scoresFile, std::ios::in);
	openck(clipFileName, clipFile, std::ios::in);
	openck(seqOutFileName,   seqOutFile, std::ios::out);
	openck(scoresOutFileName,  scoresOutFile, std::ios::out);	
	std::string line;
	std::getline(clipFile, line);

	DNASequence seq, newSeq;
	std::string title, code;
	std::string clipCode;
	ssize_t clipBegin, clipEnd;
	newSeq._ascii = 1;
	std::vector<ssize_t> qualityScores;
	std::string qualityTitle;
	while (SeqReader::GetSeq(seqFile, seq, SeqReader::noConvert)) {
		qualityScores.clear();
		ReadQualityScores(scoresFile, qualityTitle, qualityScores);
		title = seq.namestr;
		/*
			std::cout << "title: " << title << std::endl;
		*/
		ssize_t codeEndPos = title.rfind(" ");
		if (codeEndPos == title.npos) 
			codeEndPos = title.size();

		std::string code(&(title.c_str()[7]), title.c_str() + codeEndPos);
		/*
			std::cout << "searching for code: " << code << std::endl;
		*/
		clipCode = "";
		clipBegin = -1;
		clipEnd   = -1;
		while (clipFile and (clipCode == "" or clipCode < code)) {
			while (clipFile and (clipCode == "" or clipCode < code)) {
				clipFile >> clipCode >> clipBegin >> clipEnd;
				std::getline(clipFile, line);
			}
			if (!clipFile)
				break;
			if (clipCode >= code)
				break;
		}
		if (clipCode > code) {
			std::cout << "no clipping found for: " << code << std::endl;
		}
		if (clipCode == code and 
				clipBegin > 0 and clipEnd > clipBegin) {
			/*
			std::cout << "match, clipping: " << clipBegin 
								<< " " << clipEnd << std::endl;
			*/
			assert(clipBegin > 0);
			if (clipEnd > seq.length)
				clipEnd = seq.length;

			assert(clipEnd <= seq.length);
			assert(clipEnd <= qualityScores.size());
						 
			newSeq.seq = &seq.seq[clipBegin-1];
			newSeq.length = clipEnd - clipBegin + 1;
			newSeq.namestr = seq.namestr;
			newSeq.PrintSeq(seqOutFile);
			seqOutFile << std::endl;
			
			std::vector<ssize_t> clippedScores;
			clippedScores.insert(clippedScores.begin(),
													 qualityScores.begin() + clipBegin - 1,
													 qualityScores.begin() + clipEnd - 1);
			PrintQualityScores(scoresOutFile, qualityTitle, clippedScores);
		}
		else {
			std::cout << "not outputting " << code << " since clip range : " 
								<< clipBegin <<"  " << clipEnd << std::endl;
		}
	}
	return 0;
}

void ReadQualityScores(std::ifstream    &scoresIn, 
											 std::string      &title,
											 std::vector<ssize_t> &scores) {
	scoresIn >> title;
	std::string line;
	std::getline(scoresIn, line);
	ssize_t score;
	//	std::cout << "read title: " << title << std::endl;
	while (scoresIn) {
		if (scoresIn.peek() == '>') 
			break;
		scoresIn >> score;
		if (scoresIn.peek() == '\n')
			scoresIn.get();
		scores.push_back(score);
	}
	//	std::cout << "read qs of ength: " << scores.size() << std::endl;
}
		
void PrintQualityScores(std::ofstream &scoresOut,
												std::string &title,
												std::vector<ssize_t> &scores) {

	//UNUSED// ssize_t width = 30;
	ssize_t i;
	scoresOut << title << std::endl;
	for (i = 0; i < scores.size(); i++) {
		scoresOut.width(3);
		scoresOut << scores[i];
		if ((i > 0 and i % 30 == 0) or i == scores.size()-1)
			scoresOut << std::endl;
		else
			scoresOut <<" ";
	}
}

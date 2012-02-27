/***************************************************************************
 * Title:          BuildMateTable.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include <map>
//#include "regexp/boost/regex.hpp"
#include "utils.h"
#include "ParseTitle.h"
#include "MateLibrary.h"
#include "compatibility.h"
#include <regex.h>
#include "IntegralTupleStatic.h"

void PrintUsage() {
	std::cout << "usage: buildMateTable readFile ruleFile pairFile\n";
}
using namespace std;

std::map<std::string, Clone > mateMap;
int main(int argc, char* argv[]) {
	
	std::string ruleFileName;
	std::string readFileName;
	std::string pairFileName;

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}

	readFileName = argv[1];
	ruleFileName = argv[2];
	pairFileName = argv[3];
	
	std::string reportFileName = FormReportName(readFileName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	ssize_t startReadIndex = 0;

	RuleList rules;
	ParseRuleFile(ruleFileName, rules, report);
	
	std::ifstream readsIn;
	std::ofstream pairsOut;
	openck(readFileName, readsIn, std::ios::in, report);
	openck(pairFileName, pairsOut, std::ios::out, report);

	std::string line, title;
	ssize_t r;
	//	boost::smatch match;

	// support for reading multiple pair files
	ssize_t readIndex = startReadIndex;
	std::vector<std::string> mateNames;
	std::map<std::string, Clone>::iterator mapIt;						
	ssize_t numMatchedReadTitles = 0;
	ssize_t numReads = 0;
	while(std::getline(readsIn, line)) {
		if (line.size() > 1 and line[0] == '>') {
			// Found a fasta title.
			++numReads;
			if (!(ParseTitle(line, title))) {
				std::cout <<"Bad FASTA title:" << std::endl;
				std::cout << line << std::endl;
				exit(1);
			}
			ssize_t cloneFound = 0;
			//			cout << "matching " << title << endl;
			for (r = 0; r < rules.size(); r++) {
				
				_SZT_ nmatch = 3;
				regmatch_t pmatch[3];
				_INT_ mv;
				const char* titleStr = title.c_str();
				if (!(mv = regexec(&rules[r].compRegex, titleStr, nmatch, pmatch, 0))) {
					//										cout << "matched. " << r << " " <<  rules[r].type << endl;
					std::string mateBase;// = match.str(1);
						std::string matePair;// = match.str(2);
						mateBase.assign(&(titleStr[pmatch[1].rm_so]), pmatch[1].rm_eo - pmatch[1].rm_so);
						matePair.assign(&(titleStr[pmatch[2].rm_so]), pmatch[2].rm_eo - pmatch[2].rm_so);
						mateNames.push_back(mateBase);
						//						cout << " mate base: " << mateBase << endl;
						++numMatchedReadTitles;
						if ((mapIt = mateMap.find(mateBase)) == mateMap.end()) {
							Clone clone;
							clone.type = r;
							clone.ai   = readIndex;
							clone.bi   = -1;
							/*
							if (matePair == rules[r].forward)
								clone.aDir = 0;
							else
								clone.aDir = 1;
							*/
								 
							mateMap.insert(NameClonePair(mateBase, clone));
							++readIndex;
						}
						else {
							
							if ((*mapIt).second.bi >= 0) {
								std::cout << "ERROR! A clone is specified with the same mate-name three times." << std::endl;
								std::cout << (*mapIt).first << ", " << mateBase << std::endl;
								std::cout << line << std::endl;
//								exit(1);

							}
							else {
								(*mapIt).second.bi = readIndex;
								++readIndex;
								//								cout << (*mapIt).second.ai << " " << (*mapIt).second.bi << endl;
							}

							/*							if (matePair == rules[r].forward)
								(*mapIt).second.bDir = 0;
							else
								(*mapIt).second.bDir = 1;
							*/
						}
						cloneFound = 1;
						break;
						// Done searching through rules
				}
				else {
					if (mv == REG_NOMATCH) {
						//						cout << " no match " << endl;
					}
					else {
						char errorstr[1000];
						regerror(mv, &rules[r].compRegex, errorstr, 1000) ;
						cout << "error: " << errorstr << endl;
					}
				}
			} // end loop iterating over mate rules
			if (!cloneFound) {
				mateNames.push_back("");
				++readIndex;
			}
		} // end finding a read
	} // done looking thorugh all reads.
	cout << "Matched " << numMatchedReadTitles << " of " << numReads << endl;
	ssize_t maxReadIndex = readIndex;
	std::vector<Clone> clones;
	clones.resize(mateMap.size());
	cout << "writing " << maxReadIndex << " indices." << endl;
	for (readIndex = 0; readIndex < maxReadIndex; readIndex++) {
		if (mateNames[readIndex] != "") {
			mapIt = mateMap.find(mateNames[readIndex]);
			if (mapIt == mateMap.end()) {
				std::cout << "ERROR, there is an internal inconsistency.  Please" << std::endl
									<< " contact the authors." << std::endl;
				exit(1);
			}
			if ((*mapIt).second.ai != -1 and
					(*mapIt).second.bi != -1) {
				if (readIndex == (*mapIt).second.ai)
					pairsOut << (*mapIt).second.bi << " ";
				else
					pairsOut << (*mapIt).second.ai << " ";
				pairsOut << (*mapIt).second.type << std::endl;
			}
			else {
				pairsOut << "-1 -1"<< std::endl;
			}
		}
		else {
			pairsOut << "-1 -1" << std::endl;
		}
	}

	EndReport(report);
	report.close();

	return 0;
}


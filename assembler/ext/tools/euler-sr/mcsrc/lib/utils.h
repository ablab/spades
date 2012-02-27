/***************************************************************************
 * Title:          utils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _UTILS_H_
#define _UTILS_H_
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include "compatibility.h"

#define CheckedOpen


// With NOCLOBBER compiler option, files that are about to be clobbered will instead be renamed
//#define NOCLOBBER
//#undef NOCLOBBER


/* OpenIfAvailable removed: use openck with severe=-1 instead */
/*
template <typename t>
ssize_t OpenIfAvailable(std::string fileName, t &file) {
	file.open(fileName.c_str());
	if (file.good()) return 1;
	else return 0;
}
*/

// openck(fileName, file, flags, report, severe)
// report: report file ostream
//         A log entry will be written to it with the filename and flags
// severe: If the file is not openable:
//         1: error message and abort
//         0: warning message but continue.  Return 0.
//         -1: no message.  Return 0.
// return value: 1 if successfully opened, 0 otherwise
//
// Compiler variable NOCLOBBER:
// if set, files that are opened for writing will be backed up first if
// they already exist

template <typename t>
ssize_t openck(std::string fileName, t &file, 
					 std::ios_base::openmode flags= (std::ios::in),
					 std::ostream &report = std::cout,
					 ssize_t severe=1) {


#ifdef NOCLOBBER
	// In output mode (not input or append), if file already exists,
	// rename it by adding on first available suffix .v1, .v2, ..., .v99
	if ((flags & (std::ios::out | std::ios::in | std::ios::app))
			== std::ios::out) {
		struct stat buf;
		if (stat(fileName.c_str(), &buf)) {
			// File doesn't exist, so no need to rename it.
# if 0
			report << "could not stat file " << fileName << std::endl;
# endif
		} else {
			std::stringstream newFname;
			_INT_ gotname = 0;
			for (_INT_ vernum = 1; vernum < 100; vernum++) {
				newFname.str("");
				newFname << fileName << ".v" << vernum;
				if (stat(newFname.str().c_str(),&buf)) {
					gotname = 1;
					break;
				}
			}

			if (gotname && !rename(fileName.c_str(),newFname.str().c_str())) {
				report << "Renamed:         " << fileName
							 << " to " << newFname.str() << std::endl;
			} else {
				report << "Error trying to rename " << fileName;
				if (gotname) {
					report << " to " << newFname.str();
				}
				report << std::endl;
			}
		}
	}
#endif // NOCLOBBER


	// Open file and check success/failure
  file.open(fileName.c_str(), flags);
  if (!file.good()) {
		if (severe == -1) {
			return 0;
		}
    std::cout << "could not open " << fileName << std::endl;
    if (severe) 
      exit(0);
    else
     return 0;
  }

	// Log this to the report file
	if (report != file) {
		std::string openType = "Opening";
		if ((flags & (std::ios::out | std::ios::in)) == (std::ios::out | std::ios::in)) {
			openType = "Updating";
			if (flags & std::ios::app) {
				openType = "Updating (appending)";
			}
		} else if (flags & std::ios::app) {
			openType = "Appending";
		} else if (flags & std::ios::in) {
			openType = "Reading";
		} else if (flags & std::ios::out) {
			openType = "Writing";
		}

		if (flags & std::ios::ate) {
			openType += " (at end)";
		}
		if (flags & std::ios::trunc) {
			openType += " (& truncating)";
		}
		if (flags & std::ios::binary) {
			openType += " binary";
		}
		openType += ": ";
		report << std::setw(17) << std::left << openType << std::setw(0) << std::right << fileName << std::endl;
	}

  return 1;
}



template <typename T>
ssize_t Convert(std::string &inString, T &value) {
   std::stringstream instream;
   instream.str(inString);
   instream >> value;
   return (!instream.fail());
}

ssize_t PeekInput(std::ifstream &in, std::string str); 

void ParseFileName(std::string fileName, 
		   std::string &refName, 
		   std::string &qryName,
		   std::string &seqName);

//int OnFwgridNode();

class Pos {
 public:
  ssize_t rBegin, rEnd, qBegin, qEnd;
  double eValue;
};

void Tokenize(std::string &line, std::vector<std::string> &values);
ssize_t RunBlast(std::string &queryName, std::string &sbjctName, Pos &position);   
void MakeTempName(std::string &fileName, std::string ext = "");
std::string CommandLineToString(int argc, char* argv[]);
std::string CurTimeString();
void BeginReport(int argc, char *argv[], std::ofstream &out);
void EndReport(std::ofstream &out);

std::string FormReportName(std::string inName);
std::string FormReportName(char *inName);


template<typename T>
std::string NumToStr(T number) {
	std::stringstream numStrStrm;
	numStrStrm << number;
	return numStrStrm.str();
}

ssize_t IsOption(const char* arg);
ssize_t IsOption(const char* arg, const char* option);

void PrintStatus(ssize_t pos, ssize_t spacing=1000);
	
void Pause();
void WaitLock(std::string &lockFileName, _INT_ &fileDes);
void ReleaseLock(_INT_ fileDes);
#endif

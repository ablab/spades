/***************************************************************************
 * Title:          utils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "compatibility.h"
#include "utils.h"
#include <sstream>
#include <fstream>
#include <ostream>
#include <istream>
#include <fcntl.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#if defined(__sun__)
#  include <procfs.h>
#endif

using namespace std;

void ParseFileName(std::string fileName, 
		   std::string &refName, 
		   std::string &qryName,
		   std::string &seqName) {
  // Find the ending component of the path.
  ssize_t startRefInd, endRefInd, qryStartInd, qryEndInd;

  startRefInd = fileName.rfind("/");
  if (startRefInd < 0) 
    startRefInd = 0;
  
  endRefInd = fileName.find(".EN", startRefInd);
  if (endRefInd == fileName.npos) {
    std::cout << "error parsing " << fileName << std::endl;
    exit(0);
  }
  startRefInd += 1;
  refName = fileName.substr(startRefInd, endRefInd - startRefInd );

  ssize_t faStart;
  faStart = fileName.find(".fa", endRefInd);
  if (faStart == fileName.npos) {
    std::cout << "error parsing " << fileName << " " << endRefInd << std::endl;
    exit(0);
  }
  seqName = fileName.substr(endRefInd + 1, faStart  - endRefInd -1);

  qryEndInd = fileName.find(".EN", faStart);
  qryStartInd = faStart + 4;
  qryName = fileName.substr(qryStartInd, qryEndInd - qryStartInd);
}

/*
ssize_t OnFwgridNode() {
  std::string fwgridTest = "/state/partition1/test";
  std::ofstream out;
  out.open(fwgridTest.c_str());
  if (out.good()) {
    ssize_t res = system(std::string("rm "  + fwgridTest).c_str());
    return 1;
  }
  else {
    return 0;
  }
}
*/  

ssize_t RunBlast(std::string &queryName, std::string &sbjctName, Pos &position) {
  // create the temporary output name
  std::string posFile;
  MakeTempName(posFile, "locator.pos");
  std::string command;

  command = std::string(getenv("MCSRC"));
  command += "/comparative/TopHitBl2seq.pl"; 
  command += " " + queryName + " " + sbjctName + " ";
  command += posFile;

  // search the regions
  if (system(command.c_str()))
    return 0;

  // read in the resulting positions  
  std::ifstream fragPosIn;
  openck(posFile, fragPosIn, std::ios::in);
  fragPosIn >> position.rBegin 
	    >> position.rEnd
	    >> position.qBegin
	    >> position.qEnd
	    >> position.eValue;

  fragPosIn.close();
	//  int res = system(std::string("rm " + posFile).c_str());
	system(std::string("rm " + posFile).c_str());
	// TODO: interpret return status from system

  if (position.rBegin != -1 and
      position.rEnd   != -1 and
      position.qBegin != -1 and
      position.qEnd   != -1)
    return 1;
  else
    return 0;
}

void Tokenize(std::string &line, std::vector<std::string> &values) {
  std::stringstream linestrm;
  linestrm.str(line);
  std::string value;
  while (linestrm >> value) {
   values.push_back(value);
 }
}

ssize_t PeekInput(std::ifstream &in, std::string str) {
  ssize_t len = str.size();
  char *strPtr = (char*) str.c_str();
  ssize_t chkPos = 0;
  while (chkPos < len and 
	 in and  in.peek() != EOF and
	 ((char)in.peek()) == strPtr[chkPos]) {
    in.get();
   ++chkPos;
  }
  if (chkPos == len) {
    return 1;
  } 
  else {
    if (in.eof()) {
      in.clear(std::ios::eofbit);
    }
    if (!in) {
      // bad bit set, bail out
      return 0;
    }
    while (chkPos > 0) {
       in.putback(strPtr[chkPos]);
       chkPos--;
    }
    return 0;
  }  
}


void MakeTempName(std::string &fileName, std::string ext ) {
  std::string fwgridBase = "/state/partition1/";
  pid_t pid = getpid();
  std::stringstream pidstr;
  pidstr << pid;
  if (ext != "") 
    ext = "." + ext;
 
  fileName = fwgridBase + pidstr.str() + ext;
  std::ofstream testOut;
  testOut.open(fileName.c_str());
  if (testOut.good()) {
    testOut.close();
    return;
  }
  else {
    // Make in the current directory
    fileName = pidstr.str() + ext;
  }
}

std::string CommandLineToString(int argc, char* argv[]){
	ssize_t i;
	std::string result = "";
	assert(argc >= 0);
	for (i = 0; i < argc-1; i++) {
		result += argv[i];
		result += " ";
	}
	result += argv[i];
	return result;
}

std::string CurTimeString() {
	time_t curTime;
	struct tm *curTimeLocal;
	char timeBuf[100];

	if (time(&curTime) != (time_t) -1
			&& (curTimeLocal = localtime(&curTime))
			&& strftime(timeBuf, sizeof(timeBuf),
									"%a %b %e %H:%M:%S %Z %Y", curTimeLocal)) {
				return std::string(timeBuf);
	}

	return "Time not available";
}

void BeginReport(int argc, char *argv[], std::ofstream &out) {
	out << "Running:         " << argv[0] << std::endl;
	char *cwd = getenv("PWD");
	if (cwd) {
		out << "Directory:       " << cwd << std::endl;
	}
	//	time_t curTime;
	//	curTime  = time(&curTime);
	//	struct tm *curTimeTM= gmtime(&curTime);
	//	char *curTimeStr = asctime(curTimeTM);
	//	if (curTimeStr) {
	//		out << "Start:        " << curTimeStr << std::endl;
	//	}

	out << "Start:           " << CurTimeString() << std::endl;
	out << "Command:         " << CommandLineToString(argc, argv) << std::endl;
//	delete[] curTimeStr;
}

void EndReport(std::ofstream &out) {
	//	time_t curTime;
	//	curTime   = time(&curTime);
	//	struct tm *curTimeTM = gmtime(&curTime);
	//	char *curTimeStr = asctime(curTimeTM);
	//	if (curTimeStr) {
	//		out << "End:\t" << curTimeStr << std::endl;
	//	}

	out << "End:             " << CurTimeString() << std::endl;

	pid_t pid = getpid();
	ifstream in;
	stringstream statFileStrm;
	statFileStrm << "/proc/" << pid << "/status";
#if defined(__linux__)
	in.open(statFileStrm.str().c_str());
	if (in.good()) {
		_SZT_ peakBytes;
		std::string line;
		while (getline(in, line)) {
			if (sscanf(line.data(), "VmPeak:"" %lu",&peakBytes))
				break;
		}
		out << "Max Memory used: " << peakBytes << " kB" << endl;
		in.close();
	}
#elif defined(__sun__)
	in.open(statFileStrm.str().c_str(), std::ios::in | std::ios::binary);
	if (in.good()) {
		pstatus_t ps;
		in.read((char*) &ps, sizeof(pstatus_t));
		_SZT_ heap = ps.pr_brksize / 1024L;
		out << "Max Memory used: " << heap << " kB" << endl;
		in.close();
	}
#elif defined(__APPLE__)
	// TODO: OS X, and other platforms not already covered
	out << "Max Memory used: Not available" << endl;
#endif

	out << "\n-------------------------------------------------------------------------------\n\n";
}

ssize_t IsOption(const char *arg) {
	if (arg == NULL)
		return 0;
	return arg[0] == '-';
}

/******************************************************************************
 * ReportName(inName)
 *   Generate filename of report file
 *
 * If inName has format
 *       prefix + base + suffix
 *    where prefix is one of the standard output directories
 *           "" fixed/ simple/ transformed/ matetransformed/
 *    base has no slashes;
 *    and suffix is one of of the standard output file suffixes
 *       "" .contig .dot .edge .intv .iovp .mates .ovp
 *       .path .spect .sv .v .altEdges
 * then return base + "." + report
 * else return inName + "." + report
 ******************************************************************************/

std::string FormReportName(std::string inName) {
	std::string reportName = inName;

	//	string::size_type loc = reportName.find("/", 0);
	ssize_t loc = reportName.find("/", 0);
	if ( loc != reportName.npos ) {
		std::string prefix = reportName.substr(0, loc);
		if (prefix == "fixed"
				or prefix == "simple"
				or prefix == "transformed"
				or prefix == "matetransformed") {
			reportName.erase(0,loc+1);
		}
	}

	loc = reportName.rfind(".", reportName.length());
	if (loc != reportName.npos) {
		std::string suffix = reportName.substr(loc+1);
		if (suffix == "contig"
				or suffix == "dot"
				or suffix == "edge"
				or suffix == "intv"
				or suffix == "iovp"
				or suffix == "mates"
				or suffix == "ovp"
				or suffix == "path"
				or suffix == "spect"
				or suffix == "sv"
				or suffix == "v"
				or suffix == "altEdges"
				) {
			reportName.erase(loc);
		}
	}

	loc = reportName.find("/", 0);
	if (loc != reportName.npos) {
		// Filename has two slashes, don't truncate suffix/prefix
		reportName = inName;
	}

	reportName = reportName + ".report";

	return reportName;
}

std::string FormReportName(char *inName) {
	return FormReportName(std::string(inName));
}


/*****************************************************************************/

ssize_t IsOption(const char *arg, const char *option) {
	return (strcmp(arg, option) == 0);
}

void PrintStatus(ssize_t pos, ssize_t spacing) {
	ssize_t end = spacing -1;
	ssize_t line = spacing*50;
	ssize_t lineEnd = line - 1;
	if (pos % spacing == end ) 
		std::cout << "." << std::flush;
	if (pos % line == lineEnd)
		std::cout << pos + 1 << std::endl;
}

void WaitLock(std::string &lockFileName, _INT_ &lockFileDes) {
	
	struct flock fldes;
	pid_t mypid = getpid();
	fldes.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
	fldes.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
	fldes.l_start  = 0;        /* Offset from l_whence         */
	fldes.l_len    = 0;        /* length, 0 = to EOF           */
	fldes.l_pid    = mypid;
	
	FILE *f = fopen(lockFileName.c_str(), "a");
	if (f != NULL)
		fclose(f);
	else {
		std::cout << "Error creating file " << lockFileName << std::endl;
		exit(1);
	}
	std::cout << "waiting on " << lockFileName << std::endl;
	if ((lockFileDes = open(lockFileName.c_str(), O_RDWR|O_CREAT,S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH )) == -1) {
		std::cout << "error opening " << lockFileName << " " << errno << std::endl;
		exit(1);
	}
	if (fcntl(lockFileDes, F_SETLKW, &fldes) == -1) {
		std::cout << "Could not set lock on " << lockFileName << " " << errno << std::endl;
		exit(1);
	}
	std::cout << "continuing on " << lockFileName << std::endl;
}

void ReleaseLock(_INT_ lockFileDes) {
	struct flock fldes;
	pid_t mypid = getpid();
	fldes.l_type   = F_UNLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
	fldes.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
	fldes.l_start  = 0;        /* Offset from l_whence         */
	fldes.l_len    = 0;        /* length, 0 = to EOF           */
	fldes.l_pid    = mypid;
	
	if (fcntl(lockFileDes, F_WRLCK, &fldes) == -1) {
		std::cout << "could not unlock file lock" << std::endl;
		exit(0);
	}
	close(lockFileDes);
}

void Pause() {
	ssize_t a=0;
	a++;
}

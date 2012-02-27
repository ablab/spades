/***************************************************************************
 * Title:          extractseq.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <fstream>

#include <unistd.h>
#include <assert.h>
#include <sstream>
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"



void initenv(int argc, char *argv[], 
	     ssize_t & threshold,
	     ssize_t & extend,
	     ssize_t & contig_offset,
	     std::string &seqfile, 
	     std::string &positionfile,
	     std::string &outfile,
	     ssize_t &contigsPerFile);

void printusage();
void PrintSeq(std::ofstream &outfile, std::string &seq, ssize_t strlen, ssize_t lineWidth, _SSZT_ position, std::string titleContigs);

class StripLocation {
public:
  _SSZT_ start;
  _SSZT_ end;
  ssize_t  dir;
};

ssize_t ReadStrip(std::ifstream &inFile, StripLocation &loc) {
  std::stringbuf strbuf;
  //  inFile.get(strbuf);
  loc.start = -1;
  loc.end  = -1;
  if (inFile.eof())
    return 0;
  inFile >> loc.start >> loc.end;
  inFile.ignore(10000,'\n');
  if (loc.start < 0) {
    loc.dir = -1;
    loc.start = loc.start * -1;
  }
  else {
    loc.dir = 1;
  }
  if (loc.start > loc.end) {
    loc.dir = -1;
    _SSZT_ temp;
    temp= loc.start;
    loc.start = loc.end;
    loc.end   = loc.start;
  }
  return (loc.start != -1 && loc.end != -1);
}

void ReadStrips(std::ifstream &inFile, std::vector<StripLocation> &locations) {
  StripLocation loc;
  ssize_t status = 1;
  while (inFile.good() && !inFile.eof() && status != 0) {
    status = ReadStrip(inFile, loc);
    if (status) {
      locations.push_back(loc);
    }
  }
}

void ExtendStrips(std::vector<StripLocation> & locations, ssize_t extend) {
  std::vector<StripLocation>::iterator it, end;
  end = locations.end();
  for(it = locations.begin(); it != end; ++it) {
    it->end += extend;
  }
}

void GetString(std::vector<std::string> src, std::vector<_SSZT_> rl, StripLocation loc, std::string &result, std::string &titles) {
  ssize_t i;
  ssize_t curLen = 0;
  std::vector<_SSZT_> start, len;
  ssize_t startContig, endContig;
  result = "";
  std::stringstream titleStream(std::stringstream::out);
  ssize_t numContig = src.size();
  _SSZT_ genomeLength = rl[rl.size()-1];
  // rl is a vector of running lengths, the running total of the length of the sequences at this point at seq[i]
  _SSZT_ contigStartLoc, contigEndLoc;
  i = 0;

  if (loc.end > (genomeLength)) {
    std::cout << "contig not in sequence " << std::endl;
    return;
    }
  startContig = 0;
  _SSZT_ la, lb;
  while (loc.end >= rl[i] && i < numContig){
    // make sure the sequence has started before this contig
    la = rl[i+1];
    lb = rl[i];
    if (loc.start >= rl[i]){
      // start inside the contig
      startContig = i;
      contigStartLoc = loc.start - rl[i];
    }
    else {
      // start at the beginning of the contig.
      contigStartLoc = 0;  
    }

    // make sure that the sequence has started
    if (loc.start < rl[i+1] && loc.end >= rl[i]) {
      endContig = i;
      start.push_back(contigStartLoc);
    }
    
    // determine how much of this contig is extracted
    if(loc.start < rl[i+1]) {
      if (loc.end >= rl[i+1]) {
	len.push_back(src[i].size() - contigStartLoc);
      }
      else {
	la = loc.end - rl[i] - contigStartLoc;
	len.push_back(la);
      }
    }
    i++;

    if (i == rl.size() && rl[i-1] > loc.end) {
      result = "end not found";
      return;
    }
  }
    for (i = startContig; i <= endContig; i++) {
      ssize_t j;
      result.append(src[i], start[i-startContig], len[i-startContig]);
      titleStream << " contig " << i;
    }
    titles = titleStream.str();
}
int main(int argc, char *argv[]) {
  std::ifstream seqFile, locFile;
  std::ofstream outFile;
  std::string seqFileName, positionFileName, outFileName;
  std::vector<StripLocation> locations;
  ssize_t threshold = 20;
  ssize_t extend = 0;
  ssize_t contig_offset = 0;
  ssize_t contigsPerFile = -1;
  initenv(argc, argv, 
	  threshold, extend, 
	  contig_offset, 
	  seqFileName, positionFileName, outFileName,
	  contigsPerFile);
  
  _SSZT_ seqnum = 0;
  _SSZT_ startPos = 0;
  _SSZT_ genomeLength = 0;
  std::vector<_SSZT_> runningLength, revRunningLength;
  std::vector<std::string>sequences, revSequences;
  DNASequence *seq, *rev;
  //  char *seq, *rev, *seqname;
  _SSZT_ len, s;

  locFile.open(positionFileName.c_str());
  ReadStrips(locFile, locations);
  ExtendStrips(locations, extend);

  seqFile.open(seqFileName.c_str());
  SeqReader::MaskRepeats();
  SeqReader reader(&seqFile);
  startPos = 0;
  reader.GetSeq(seq);
  seqFile.close();
  rev = new DNASequence;
  MakeRC(*seq, *rev);


  DNASequence subSeq;
  outFile.open(outFileName.c_str());
  std::vector<StripLocation>::iterator it, end;
  for (it = locations.begin(); it != locations.end(); ++it) {
    if (it->dir == 1) 
      subSeq.seq = &seq->seq[it->start];
    else 
      subSeq.seq = rev->seq;
    subSeq.length = it->end - it->start + 1;
    std::stringstream namestr;
    namestr << outFileName << "_" << it->start << "_" << it->end;
    std::cout << "subseq length: " << subSeq.length << std::endl;
    subSeq.StoreName((char*)namestr.str().c_str());
    subSeq.PrintSeq(outFile);
  }

  outFile.close();
  return 0;
}


void PrintSeq(std::ofstream &outfile, std::string &seq, ssize_t strlen, ssize_t lineWidth, _SSZT_ position, std::string titleContigs) {
  ssize_t lenPrinted = 0;
  ssize_t i, end;
  char *substr =(char*) new char[lineWidth+1];
  outfile  << ">" << position << " " << titleContigs << std::endl;
  while (lenPrinted < strlen) {
    end = (lenPrinted+lineWidth > strlen) ? strlen : lenPrinted + lineWidth;
    strncpy(substr, seq.c_str() + lenPrinted, end - lenPrinted);
    substr[end-lenPrinted] = '\0';
    //    std::cout << "got : " << substr << std::endl;
    outfile << substr << std::endl;
    lenPrinted += lineWidth;
  }
}

void initenv(int argc, char *argv[], ssize_t & threshold, ssize_t & extend, ssize_t & contig_offset, std::string &inpfile, std::string &locFile,std::string &outfile,
	     ssize_t &contigsPerFile) {
  ssize_t copt;
  inpfile = "?";
  locFile = "?";
  outfile = "?";
  while ( (copt=getopt(argc, argv, "i:o:l:e:p:c")) != EOF) {
    switch(copt) {
    case 'c':
      contigsPerFile = 1;
      continue;
    case 'i':
      inpfile = optarg;
      continue;
    case 'o':
      outfile = optarg;
      continue;
    case 'l':
      locFile= optarg;
      continue;
    case 't':
      threshold = atoi(optarg);
      continue;
    case 'e':
      extend = atoi(optarg);
      continue;
    case 'p':
      contig_offset = atoi(optarg);
      continue;
    default:
      printusage();
      exit(1);
    }
  }
  if (inpfile == "?"|| locFile == "?" || outfile == "?") {
    std::cout << "inp " << inpfile << " out " << outfile << " loc " << locFile << std::endl;
    printusage();
    exit(1);
  }
}
void printusage() {
  std::cout << "usage:  mergeblocks  -i seq1  [-d distThresh] -o output [-t] " << std::endl;
  std::cout << "   -i input       : sequence file." << std::endl;
  std::cout << "   -l locations   : locations file, format: start end dir " << std::endl;
  std::cout << "   -o output      : output file " << std::endl;
  std::cout << "   -e extend      : amount to extend each strip (add on) " << std::endl;
  std::cout << "   -p probelen    : amount to adjust each read length by (probe-length) " << std::endl;
}

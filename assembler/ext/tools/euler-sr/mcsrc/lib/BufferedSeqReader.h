/***************************************************************************
 * Title:          BufferedSeqReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/07/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "utils.h"
#include "SeqReader.h"
#include <fstream>
#include <vector>
using namespace std;
template<ssize_t BUF_SIZE>
class BufferedSeqReader {
public:
  std::ifstream *seqIn;
	std::vector<DNASequence> seqBuffer;
  ssize_t cur;
  ssize_t end;
	FILE* lockfile;
	string lockFileName;
	ssize_t lockFileDes;
	void Init() {
    Recharge();
		lockFileName = "";
		lockFileDes = -1;
	}		
  void Init(std::string fileName) {
		seqIn = new std::ifstream;
    seqBuffer.resize(BUF_SIZE);
		Init();
  } 
	void Init(std::ifstream *in) {
		seqIn = in;
		Init();
	}

	ssize_t Recharge() {
		return Recharge(seqBuffer);
	}

  ssize_t Recharge(std::vector<DNASequence> &newSeqBuffer) {
    ssize_t curRead = 0;
		DNASequence read;
		if (newSeqBuffer.size() < BUF_SIZE) 
			newSeqBuffer.resize(BUF_SIZE);

    while (curRead < BUF_SIZE and SeqReader::GetSeq(*seqIn, read, SeqReader::noConvert) ){
			newSeqBuffer[curRead] = read;
			assert(curRead < newSeqBuffer.size());
			curRead++;
		}
    end = curRead;
    cur = 0;
    return end;
  }

  ssize_t GetSeq(DNASequence &seq) {
   if (cur == end and Recharge() == 0)
     return 0;
   seq = seqBuffer[cur];
   cur++;
   return 1;
  }
	void Reset() {
		seqIn->close();
		seqIn->clear();
	}


	ssize_t GetRead(DNASequence &seq) {
		if (cur == end and RechargeRead() == 0) 
			return 0;
		seq = seqBuffer[cur];
		cur++;
		return 1;
	}
  ssize_t RechargeRead() {
    ssize_t curRead = 0;
		DNASequence read;
    while (curRead < BUF_SIZE and SeqReader::GetRead(*seqIn, read) ){
			seqBuffer[curRead] = read;
			assert(curRead < seqBuffer.size());
			curRead++;
		}
    end = curRead;
    cur = 0;
    return end;
  }	
};

  

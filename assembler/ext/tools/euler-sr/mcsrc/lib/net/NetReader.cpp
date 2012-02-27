/***************************************************************************
 * Title:          NetReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "NetReader.h"
#include "NetFile.h"


#include "utils.h"

ssize_t NetReader::ParseLevel(std::ifstream &in, ssize_t &level) {
  level = 0;
  while(in && in.peek() == ' ') {
    in.get();
    level++;
  }
  return in.good();
}

void NetReader::ParseKeywordValues(std::ifstream &in, Net &net) {

  std::string key, type;

  while (in && in.peek() != '\n') {
    in >> key;
    if (key == "id") {
      in >> net.id;
    }
    else if (key ==  "score") {
      in >> net.score;
    }
    else if (key ==  "ali") {
      in >> net.ali;
    }
    else if (key ==  "qFar") {
      in >> net.qFar;
    }
    else if (key ==  "qOver") {
      in >> net.qOver;
    }
    else if (key ==  "qDup") {
      in >> net.qDup;
    }
    else if (key ==  "type") {
      in >> type;
      if (type == "top") {
	net.type = Net::top;
      }
      else if (type == "syn") {
	net.type = Net::syn;
      }
      else if (type == "inv") {
	net.type = Net::inv;
      }
      else if (type == "nonSyn") {
	net.type = Net::nonSyn;
      }
    }
    else if (key ==  "tN") {
      in >> net.tN;
    }
    else if (key ==  "qN") {
      in >> net.qN;
    }
    else if (key ==  "tR") {
      in >> net.tR;
    }
    else if (key ==  "qR") {
      in >> net.qR;
    }
    else if (key ==  "tNewR") {
      in >> net.tNewR;
    }
    else if (key ==  "qNewR") {
      in >> net.qNewR;
    }
    else if (key ==  "tOldR") {
      in >> net.tOldR;
    }
    else if (key ==  "qOldR") {
      in >> net.qOldR;
    }
    else if (key ==  "tTrf") {
      in >> net.tTrf;
    }
    else if (key ==  "qTrf") {
      in >> net.qTrf;
    }
    /*
          o top -- Chain is top-level, not a gap filler
          o syn -- Chain is on same chromosome and in same direction as parent
          o inv -- Chain is on same chromosome on opposite direction from parent
          o nonSyn -- Chain is on a different chromosome from parent 
    */

  }


}
void NetReader::ReadNetFile(std::string &inFileName, NetFile &netFile) {
  std::ifstream file;
  openck(inFileName, file);

  std::string tempstr;
  file >> tempstr;
  file >> netFile.name;
  file >> netFile.length;
  file.get(); // discard the newline
  ssize_t level;
  Net* net;
  std::string cls, orientation;
  while(file) {
    if (!ParseLevel(file, level))
      break;
    
    net = new Net;
    net->level = level;
    if (!(file >> cls >> net->tStart >> net->tSize 
	  >> net->chrom >> orientation >> net->qStart >> net->qSize)){
      delete net;
      break;
    }
    if (cls == "fill") 
	net->chainClass = Net::fill;
    else if (cls == "gap") 
	net->chainClass = Net::gap;

    if (orientation == "+") 
      net->orientation = 0;
    else
      net->orientation = 1;

    ParseKeywordValues(file, *net);
    file.get(); // discard the newline/
    netFile.nets.push_back(net);
   }

}

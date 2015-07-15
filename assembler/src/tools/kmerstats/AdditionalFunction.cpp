//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "AdditionalFunction.h"

#include <map>
#include <fstream>
#include <iostream>

void createDatFile(std::map<char, int> &e) {
    std::ofstream out("daily.dat");
    if(!out) {
        std::cout << "Cannot open file.\n";
        return;
    }
    std::map <int, int> dt;

    for (std::map<char,int>::iterator it = e.begin() ; it != e.end(); ++it) {
        if (dt.find((*it).second) == dt.end()) {
            dt[(*it).second] = 1;
        } else {
            ++dt[(*it).second];
        }
    }

    for (std::map<int,int>::iterator it = dt.begin() ; it != dt.end(); ++it) {
        out << (*it).first << " " << (*it).second << "\n";
    }
}

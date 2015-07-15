//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <string>
#include <fstream>

using namespace std;

namespace dipspades {

vector<string> GetAllLinesFromFile(string filename);

string cut_fname_from_path(string path);

bool fname_valid(string fname);

}

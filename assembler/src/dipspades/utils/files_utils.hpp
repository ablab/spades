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

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "files_utils.hpp"

namespace dipspades {

vector<string> GetAllLinesFromFile(string filename){
    ifstream freader(filename.c_str());
    vector<string> lines;
    if(!freader.fail()){
        while(!freader.eof()){
            string new_line;
            getline(freader, new_line);
            if(new_line != "")
                lines.push_back(new_line);
        }
    }
    return lines;
}

string cut_fname_from_path(string path){
    string res;
    for(size_t i = path.size() - 1; i > 0; i--)
        if(path[i] == '/'){
            res = path.substr(i + 1, path.size() - i - 1);
            break;
        }

    for(size_t i = res.size() - 1; i > 0 ; i--)
        if(res[i] == '.'){
            res = res.substr(0, i);
            break;
        }

    return res;
}

bool fname_valid(string fname){
    ifstream out(fname.c_str());
    return !out.fail();
}

}

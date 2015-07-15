//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "standard_vis.hpp"

namespace online_visualization {

    bool IsNumber(const string& s) {
         if (s.empty())
             return false;
         for  (auto iter = s.begin(); iter != s.end(); ++iter) {
            if (!std::isdigit(*iter))
                return false;
         }
         return true;
    }
            
    int GetInt(string str) {
        stringstream ss(str);
        int ans;
        ss >> ans;
        return ans;
    }
    
    vector<string> SplitInTokens(stringstream& args) { 
        vector<string> answer;
        while (!args.eof()) {
            string arg;
            args >> arg;
            answer.push_back(arg);
        }
        return answer;
    }
}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <string>
#include <memory.h>

using namespace std;

namespace dipspades {

template<class T>
class LCSCalculator{

    int ** mask;

    void Initialize(size_t length1, size_t length2){
        mask = new int*[length1 + 1];
        for(size_t i = 0; i < length1 + 1; i++){
            mask[i] = new int[length2 + 1];
            memset(mask[i], 0, sizeof(int) * (length2 + 1));
        }
    }

    void Terminate(size_t length){
        for(size_t i = 0; i < length + 1; i++)
            delete[] mask[i];
        delete[] mask;
    }

    void RecursiveLCSLenCalc(vector<T> str1, vector<T> str2, size_t i1, size_t i2){

//        cout << i1 << " - " << i2 << "; " << str1[i1 - 1] << " - " << str2[i2 - 1] << "; ";

        if(str1[i1 - 1] == str2[i2 - 1]){
//            cout << "1st case; ";
            mask[i1][i2] = mask[i1 - 1][i2 - 1] + 1;
        }
        else{

//            cout << "2nd case; ";
            int res1 = mask[i1][i2 - 1];
            int res2 = mask[i1 - 1][i2];

            mask[i1][i2] = max<int>(res1, res2);
        }
    }

    int LCSLengthCalculation(vector<T> str1, vector<T> str2){

        for(size_t i = 1; i <= str1.size(); i++)
            for(size_t j = 1; j <= str2.size(); j++){
                RecursiveLCSLenCalc(str1, str2, i, j);
            }

        return mask[str1.size()][str2.size()];
    }

    vector<T> RecursiveRestoreLCS(vector<T> str1, vector<T> str2,
            size_t i, size_t j){
        vector<T> res;
        if(i == 0 || j == 0){
            return res;
        }

        if(str1[i - 1] == str2[j - 1]){
            res = RecursiveRestoreLCS(str1, str2, i - 1, j - 1);
            res.push_back(str1[i - 1]);
            return res;
        }

        if(mask[i][j - 1] > mask[i - 1][j])
            return RecursiveRestoreLCS(str1, str2, i, j - 1);
        else
            return RecursiveRestoreLCS(str1, str2, i - 1, j);
    }

    vector<T> RestoreLCS(vector<T> string1, vector<T> string2, size_t){

//        cout << "LCS string length - " << lcs_length << endl;
        vector<T> lcs = RecursiveRestoreLCS(string1, string2, string1.size(), string2.size());
        return lcs;
    }

public:

    vector<T> LCS(vector<T> string1, vector<T> string2){
        vector<T> res;
        if(string1.size() == 0 || string2.size() == 0)
            return res;

        Initialize(string1.size(), string2.size());

        int lcs_length = LCSLengthCalculation(string1, string2);
        res = RestoreLCS(string1, string2, lcs_length);
        Terminate(string1.size());

        return res;
    }

    vector<size_t> GetPosVectorFromLeft(vector<T> string, vector<T> lcs){
        vector<size_t> pos;

        if(string.size() == 0 || lcs.size() == 0)
            return pos;

        int str_ind = 0;
        for(size_t i = 0; i < lcs.size(); i++){
            while(string[str_ind] != lcs[i]){
                str_ind++;
            }
            pos.push_back(str_ind);
            str_ind++;
        }

        VERIFY(lcs.size() == pos.size());

        return pos;
    }

    vector<size_t> GetPosVector(vector<T> string, vector<T> lcs){
        vector<size_t> pos;
        if(string.size() == 0 || lcs.size() == 0)
            return pos;

        int lcs_ind = int(lcs.size() - 1);
        int str_size = int(string.size());
        for(int i = str_size - 1; i >= 0 && lcs_ind >= 0; i--)
            if(string[i] == lcs[lcs_ind]){
                pos.insert(pos.begin(), size_t(i));
                lcs_ind--;
            }

        VERIFY(pos.size() == lcs.size());

        return pos;
    }
};

}

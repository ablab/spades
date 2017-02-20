//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <vector>
#include "utils/simple_tools.hpp"

/*
 * Little modified copy-paste from http://www.merriampark.com/ldcpp.htm
 */
inline size_t edit_distance(const std::string &source, const std::string &target) {

    // Step 1

    const size_t n = source.length();
    const size_t m = target.length();
    if (n == 0) {
        return m;
    }
    if (m == 0) {
        return n;
    }

    // Good form to declare a TYPEDEF

    typedef std::vector<std::vector<size_t> > Tmatrix;

    Tmatrix matrix(n + 1);

    // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
    // allow for allocation on declaration of 2.nd dimension of vec of vec

    for (size_t i = 0; i <= n; i++) {
        matrix[i].resize(m + 1);
    }

    // Step 2

    for (size_t i = 0; i <= n; i++) {
        matrix[i][0] = i;
    }

    for (size_t j = 0; j <= m; j++) {
        matrix[0][j] = j;
    }

    // Step 3

    for (size_t i = 1; i <= n; i++) {

        const char s_i = source[i - 1];

        // Step 4

        for (size_t j = 1; j <= m; j++) {

            const char t_j = target[j - 1];

            // Step 5

            size_t cost;
            if (s_i == t_j) {
                cost = 0;
            }
            else {
                cost = 1;
            }

            // Step 6

            const size_t above = matrix[i - 1][j];
            const size_t left = matrix[i][j - 1];
            const size_t diag = matrix[i - 1][j - 1];
            size_t cell = std::min(above + 1, std::min(left + 1, diag + cost));

            // Step 6A: Cover transposition, in addition to deletion,
            // insertion and substitution. This step is taken from:
            // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
            // Enhanced Dynamic Programming ASM Algorithm"
            // (http://www.acm.org/~hlb/publications/asm/asm.html)

            if (i > 2 && j > 2) {
                size_t trans = matrix[i - 2][j - 2] + 1;
                if (source[i - 2] != t_j) trans++;
                if (s_i != target[j - 2]) trans++;
                if (cell > trans) cell = trans;
            }

            matrix[i][j] = cell;
        }
    }

    // Step 7

    return matrix[n][m];
}

inline std::pair<std::pair<int, int>, std::string> best_edit_distance_cigar(const std::string &source,
                                                                            const std::string &target) {

    // Step 1

    const size_t n = source.length();
    const size_t m = target.length();
//  if (n == 0) {
//    return m;
//  }
//  if (m == 0) {
//    return n;
//  }

    // Good form to declare a TYPEDEF

    typedef std::vector<std::vector<int> > Tmatrix;

    Tmatrix matrix(n + 1);

    // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
    // allow for allocation on declaration of 2.nd dimension of vec of vec

    for (size_t i = 0; i <= n; i++) {
        matrix[i].resize(m + 1);
    }

    // Step 2

    for (size_t i = 0; i <= n; i++) {
        matrix[i][0] = (int) i;
    }

    for (size_t j = 0; j <= m; j++) {
        matrix[0][j] = 0; //free inserts in front
    }

    // Step 3

    for (size_t i = 1; i <= n; i++) {

        const char s_i = source[i - 1];

        // Step 4

        for (size_t j = 1; j <= m; j++) {

            const char t_j = target[j - 1];

            // Step 5

            int cost;
            if (s_i == t_j) {
                cost = 0;
            }
            else {
                cost = 1;
            }

            // Step 6

            const int above = matrix[i - 1][j];
            const int left = matrix[i][j - 1];
            const int diag = matrix[i - 1][j - 1];
            int cell = std::min(above + 1, std::min(left + 1, diag + cost));

            // Step 6A: Cover transposition, in addition to deletion,
            // insertion and substitution. This step is taken from:
            // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
            // Enhanced Dynamic Programming ASM Algorithm"
            // (http://www.acm.org/~hlb/publications/asm/asm.html)

//      if (i>2 && j>2) {
//        int trans=matrix[i-2][j-2]+1;
//        if (source[i-2]!=t_j) trans++;
//        if (s_i!=target[j-2]) trans++;
//        if (cell>trans) cell=trans;
//      }

            matrix[i][j] = cell;
        }
    }

    // Step 7
    int min = matrix[n][m];
    size_t min_m = m;

    for (size_t j = 0; j <= m; j++) {
        if (min > matrix[n][j]) {
            min = matrix[n][j];
            min_m = j;
        }
    }

//  INFO("min = "<<min<< " min_m = "<< min_m);
    std::string res = "";
    char last_operation = 0;
    int cnt_last_operation = 0;
    size_t cur_pos_i = n;
    size_t cur_pos_j = min_m;
    char cur_operation = 0;


//  if (min > 0) {
//      for (int i = 0; i <= n; i++) {
//        INFO(ToString(matrix[i]));
//      }
//  }

    while ((cur_pos_i > 0) && (cur_pos_j > 0)) {
        if (matrix[cur_pos_i - 1][cur_pos_j] < matrix[cur_pos_i][cur_pos_j]) {
            cur_operation = 'I';
            cur_pos_i--;
        }
        else {
            if (matrix[cur_pos_i][cur_pos_j - 1] < matrix[cur_pos_i][cur_pos_j]) {
                cur_operation = 'D';
                cur_pos_j--;
            }
            else {
                cur_operation = 'M';
                cur_pos_i--;
                cur_pos_j--;
            }
        }
        if (cur_operation != last_operation) {
            if (last_operation != 0)
                res = ToString(cnt_last_operation) + last_operation + res;
            last_operation = cur_operation;
            cnt_last_operation = 1;
        }
        else {
            cnt_last_operation++;
        }
    }
    res = ToString(cnt_last_operation) + last_operation + res;
    return std::make_pair(std::make_pair(cur_pos_j, min_m), res);
}

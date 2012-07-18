//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
  * union.hpp
  *
  *  Created on: 11.05.2011
  *      Author: snikolenko
  */

#ifndef CPCOUNT_UNION_H
#define CPCOUNT_UNION_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstdarg>
#include <algorithm>

using namespace std;

inline bool zero_size(const vector<int> & v) {
  return v.size() == 0;
}

//Straight out of Cormen
class unionFindClass {
public:
  unionFindClass(size_t size) : data(size, -1), rank(size, 0), sizes(size, 0) {}

  void unionn(int x, int y) {
    link(find_set(x),find_set(y));
  }

  int find_set(int x) const {
    if (data.at(x) == -1) {
      data.at(x) = x;
      sizes.at(x) = 1;
    } else if (data.at(x) != x) {
      data[x] = find_set(data[x]);
    }

    return data.at(x);
  }

  size_t num_classes() {
    size_t count = 0;
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i] == (int)i) count++;
    }
    return count; //return *max_element(data.begin(), data.end());
  }

  int num_elements() {
    return count_if(data.begin(), data.end(), bind2nd(not_equal_to<int>(), -1));
  }

  int size() {
    return data.size();
  }

  int operator[](size_t i) {
    return data.at(i);
  }

  size_t set_size(int i) const {
    if (data.at(i) == -1)
      return 0;

    int el = find_set(i);
    return sizes.at(el);
  }
  
  void get_classes (vector<vector<int > > & otherWay) const {
    otherWay.resize(data.size());
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i] != -1) 
        otherWay.at(find_set(data[i])).push_back(i);
    }
    otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), zero_size), otherWay.end());
  }

  void get_classes (vector<int> & consolidatedData, vector<vector<int > > & otherWay) {
    otherWay.resize(data.size());
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i] != -1) 
        otherWay.at(find_set(data[i])).push_back(i);
    }
    otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), zero_size), otherWay.end());

    consolidatedData.resize(data.size(), -1);
    for (size_t i = 0; i < otherWay.size(); i++)
      for (size_t j = 0; j < otherWay[i].size(); j++)
        consolidatedData.at(otherWay[i][j]) = i;
  }

private:
  void link (int x, int y) {
    if (rank.at(x) > rank.at(y)) {
      data.at(y) = x;
    } else {
      data.at(x) = y;
      if (rank.at(x) == rank.at(y))
        rank.at(y)++;
    }
    unsigned sum = sizes[x] + sizes[y];
    sizes[x] = sizes[y] = sum;
  }

  mutable std::vector<int> data;
  std::vector<int> rank;
  mutable std::vector<unsigned> sizes;
};

#endif

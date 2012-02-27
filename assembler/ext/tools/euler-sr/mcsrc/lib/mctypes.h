/***************************************************************************
 * Title:          mctypes.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MCTYPES_H_
#define MCTYPES_H_

#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <iostream>

typedef std::vector<std::string> StringVector;
typedef std::set<std::string> StringSet;
typedef std::vector<ssize_t> BitVector;
typedef std::vector<double> FloatVector;
typedef std::vector<FloatVector> FloatMatrix;

typedef std::vector<ssize_t> IntVector;
typedef std::vector<IntVector> IntMatrix;



template <typename t>
void ClearMatrix(std::vector<std::vector< t> > &matrix) {
  ssize_t i;
  for (i = 0; i < matrix.size(); i++ ) 
    matrix[i].clear();
  matrix.clear();
}

template <typename t>
void CreateMatrix(std::vector<std::vector< t> > &matrix, ssize_t rows, ssize_t cols) {
  ssize_t r, c;
  matrix.resize(rows);
  for (r = 0; r < rows; r++) {
    matrix[r].resize(cols);
    for (c = 0; c < cols; c++) 
      matrix[r][c] = 0;
  }
}

template <typename t>
void PrintMatrix(std::vector<std::vector< t> > &matrix, std::ostream &out, ssize_t width=5) {
  ssize_t r, c;
  for (r = 0; r < matrix.size(); r++) {
    for (c = 0; c < matrix[r].size() - 1; c++) {
      out.width(width);
      out << matrix[r][c] << " ";
    }
    if (matrix[r].size() > 0) {
      out.width(width);
      out << matrix[r][c] << std::endl;
    }
  }
}

template <typename t>
void ReadMatrix(std::istream &in, std::vector<std::vector< t> > &matrix) {
  std::vector< t > row;
  std::string line;
  t value;

  while (std::getline(in, line)) {
    std::stringstream strin;
    strin.str(line);
    row.clear();
    while ( strin >> value ) {
      row.push_back(value);
    }
    matrix.resize(matrix.size() + 1);
    matrix[matrix.size()-1] = row;
  }
}


    
    
    
	
    
      
  

#endif

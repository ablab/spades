/*
 * union.h
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#ifndef CPCOUNT_UNION_H
#define CPCOUNT_UNION_H

#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;
bool pred1(const vector<int> & v) {return v.size() == 0;}

//Straight out of Cormen
class unionFindClass {
public:
	unionFindClass(int size) : data(size, -1), rank(size, 0) {}

	void unionn( const int & x, const int & y) {
		link(find_set(x),find_set(y));
	}

	int find_set(const int & x) {
		if (data.at(x) == -1) {
			data.at(x) = x;
		} else if (data.at(x) != x) {
			data[x] = find_set(data[x]);
		}
		return data.at(x);
	}

	int num_classes() {
		int count = 0;
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

	int operator[](int i) {
		return data.at(i);
	}

	void get_classes (vector<vector<int > > & otherWay) {
		otherWay.resize(data.size());
		for (size_t i = 0; i < data.size(); i++) {
			if (data[i] != -1) 
				otherWay.at(find_set(data[i])).push_back(i);
		}
		otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), pred1), otherWay.end());
	}

	void get_classes (vector<int> & consolidatedData, vector<vector<int > > & otherWay) {
		otherWay.resize(data.size());
		for (size_t i = 0; i < data.size(); i++) {
			if (data[i] != -1) 
				otherWay.at(find_set(data[i])).push_back(i);
		}
		otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), pred1), otherWay.end());

		consolidatedData.resize(data.size(), -1);
		for (size_t i = 0; i < otherWay.size(); i++) 
			for (size_t j = 0; j < otherWay[i].size(); j++) 
				consolidatedData.at(otherWay[i][j]) = i;		
	}



private:
	void link (const int & x, const int & y) {
		if (rank.at(x) > rank.at(y)) {
			data.at(y) = x;
		} else {
			data.at(x) = y;
			if (rank.at(x) == rank.at(y)) rank.at(y)++;
		}
	}

	vector<int> data;
	vector<int> rank;

};


#endif 

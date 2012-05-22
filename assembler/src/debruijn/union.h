#ifndef CPCOUNT_UNION_H
#define CPCOUNT_UNION_H

#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <list>
#include <queue>
#include <cstdarg>
#include <algorithm>

bool pred1(const vector<int> & v) {return v.size() == 0;}
//TODO this class seriously needs to be refactored . 
class UnionFindClass {
    public:
        UnionFindClass(int size) : data(size, -1), rank(size, 0) {
            for(int i = 0 ; i < size ; i++)
                unionn(i,i);
        }

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
            for (int i = 0; i < (int)data.size(); i++) {
                if (data[i] == i) count++;
            }
            return count; 
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

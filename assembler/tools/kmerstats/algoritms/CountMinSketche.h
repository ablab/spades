#ifndef COUNTMINSKETCHE_H
#define COUNTMINSKETCHE_H

#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <algorithm>

int hashFunction(char i, int w) {
    return i % w;
}


void CountMinSketch(char *p, off_t size, double eps, double sigma) {
    int n = 0;
    int delta = 0;
    //std::map<char, int> m_count;
    int w = sigma / eps + 1;
    int d = log(1/sigma) + 1;

    int **m_count = new int *[d];
    for (int i = 0; i < d; i++){
        m_count[i] = new int [w];
    }
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < w; ++j) {
            m_count[i][j] = 0;
        }
    }


    for (off_t len = 0; len < size; ++len) {
        ++n;
        for (int j = 0; j < d; j++) {
            std::cout << "sss " << hashFunction(p[len], w) << std::endl;
            m_count[j, hashFunction(p[len], w)];
        }
    }

    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < w; ++j) {
     //       m_count[i][j] = 0;
           // std::cout << m_count[i][j] << std::endl;
        }
    }
/*
    std::map<char,int>::iterator it;
    for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
        std::cout << (*it).first << " => " << (*it).second << std::endl;

    }*/

    for (int i = 0; i < d; i++) {
        delete[] m_count[i];
    }
    delete[] m_count;
}


#endif // COUNTMINSKETCHE_H

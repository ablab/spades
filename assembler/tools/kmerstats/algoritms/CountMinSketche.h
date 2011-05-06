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
void CountMinSketch(char *p, off_t size, double eps, double sigma) {
    int w = 2.718281/eps + 1;
    int d = log(1 / sigma) + 1;

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
        for (int j = 0; j < d; j++) {
            int r = p[len] % w;
            m_count[j][r] += 1;
        }
    }

    int t = 's' % w;
    int min = m_count[0][t];
    for (int i = 0; i < d; ++i) {
        if (min < m_count[i][t]) {
            min = m_count[i][t];
        }
    }

    std::cout << min << std::endl;

    for (int i = 0; i < d; i++) {
        delete[] m_count[i];
    }
    delete[] m_count;
}


#endif // COUNTMINSKETCHE_H

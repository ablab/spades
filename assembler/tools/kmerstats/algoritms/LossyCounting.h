#ifndef LOSSYCOUNTING_H
#define LOSSYCOUNTING_H

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


void LossyCounting(char *p, off_t size, int k) {
    int n = 0;
    int delta = 0;
    std::map<char, int> m_count;

    for (off_t len = 0; len < size; ++len) {
        n++;
        if (m_count.find(p[len]) == m_items.end()) {
            ++m_count[ p[len] ];
        } else {
            m_count[ p[len] ] = 1 + delta;
        }

        if (n/k != delta) {
            delta = n/k;

            std::map<char,int>::iterator it;
            for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
                if ((*it).second < delta) {
                    m_count.erase(it);
                }
            }
        }
    }

    std::cout << m_count.size() << std::endl;
    std::map<char,int>::iterator it;
    for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
        std::cout << (*it).first << " => " << (*it).second << std::endl;

    }
}

#endif // LOSSYCOUNTING_H

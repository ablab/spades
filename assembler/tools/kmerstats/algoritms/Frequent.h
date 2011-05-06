#ifndef FREQUENT_H
#define FREQUENT_H

#include <iostream>
#include <map>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <algorithm>

void Frequent(char *p, off_t size, int k) {
    int n = 0;
    std::map<char, int> m_count;
    for (off_t len = 0; len < size; ++len) {
        n++;
        if (m_count.find(p[len]) == m_items.end()) {
            ++m_count[ p[len] ];
        } else if (m_items.size() < k - 1) {
            m_count[ p[len] ] = 1;
        } else {
            std::map<char,int>::iterator it;
            for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
                --(*it).second;
                if ((*it).second == 0) {
                    m_count.erase(it);
                }
            }
        }
    }

    std::map<char,int>::iterator it;
    for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
        std::cout << (*it).first << " => " << (*it).second << std::endl;

    }
}


#endif // FREQUENT_H

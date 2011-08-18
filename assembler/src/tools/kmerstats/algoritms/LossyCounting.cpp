#include "LossyCounting.h"
#include "../AdditionalFunction.h"

void WrapperLossyCounting(char *p, off_t size, int kmer, double k, double l) {
    LossyCounting(p, size, 10);
}


void LossyCounting(char *p, off_t size, int k) {
    int n = 0;
    int delta = 0;
    std::map<char, int> m_count;

    for (off_t len = 0; len < size; ++len) {
        n++;
        if (m_count.find(p[len]) == m_count.end()) {
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

    createDatFile(m_count);
}



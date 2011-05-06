#include "SpaceSaving.h"
typedef std::pair<char, int> PairType;
struct CompareSecond
{
    bool operator()(const PairType& left, const PairType& right) const
    {
        return left.second < right.second;
    }
};

void SpaceSaving(char *p, off_t size, int k) {
    int n = 0;
    std::vector<char> m_items;
    std::map<char, int> m_count;

    for (off_t len = 0; len < size; ++len) {
        n++;
        if (std::find(m_items.begin(), m_items.end(), p[len]) == m_items.end()) {
            ++m_count[ p[len] ];
        } else if (m_items.size() < k ) {
            m_items.push_back(p[len]);
            m_count[ p[len] ] = 1;
        } else {
            std::pair<char, int> min = *std::min_element(m_count.begin(), m_count.end(), CompareSecond());
            m_count[ p[len] ] = min.second + 1;
            m_count.erase(m_count.find(min.first));
        }
    }

    std::cout << m_count.size() << std::endl;
    std::map<char,int>::iterator it;
    for ( it = m_count.begin() ; it != m_count.end(); ++it ) {
        std::cout << (*it).first << " => " << (*it).second << std::endl;

    }
}

void WrapperSpaceSaving(char *p, off_t size, int kmer, double k = 0.0, double l = 0.0) {
    SpaceSaving(p, size, 30);
}

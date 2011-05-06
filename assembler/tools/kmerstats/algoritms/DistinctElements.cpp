#include "DistinctElements.h"
#include <cmath>

void WrapperDistinctElements(char *p, off_t size, double k, double l) {
    DistinctElements(p, size);
}

int zeros(int p) {
    int count = 0;
    int x;
    while (1) {
        x = p % 2;
        if (x != 0 || p == 0)
            break;
        p = p / 2;
        ++count;
    }
    return count;
}

void DistinctElements(char *p, off_t size) {
    int z = 0;
    for (off_t len = 0; len < size; ++len) {
        int r = p[len] % size;
        if (zeros( r ) > z) {
            z = zeros( r );
        }
    }

    std::cout << pow(2, z + 0.5) << std::endl;
}


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

void WrapperLossyCounting(char *p, off_t size, int kmer, double k = 0.0, double l = 0.0);

void LossyCounting(char *p, off_t size, int k);

#endif // LOSSYCOUNTING_H

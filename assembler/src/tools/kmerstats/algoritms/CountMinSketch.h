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


void WrapperCountMinSketch(char *p, off_t size, int kmer, double k, double l);

void CountMinSketch(char *p, off_t size, double eps, double sigma, int kmer);

#endif // COUNTMINSKETCHE_H

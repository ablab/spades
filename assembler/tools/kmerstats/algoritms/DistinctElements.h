#ifndef DISTINCTELEMENTS_H
#define DISTINCTELEMENTS_H


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


void WrapperDistinctElements(char *p, off_t size, double k, double l);

void DistinctElements(char *p, off_t size);

#endif // DISTINCTELEMENTS_H

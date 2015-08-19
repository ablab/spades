//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

void WrapperFrequent(char *p, off_t size, int kmer, double k, double l);

void Frequent(char *p, off_t size, int k);

#endif // FREQUENT_H

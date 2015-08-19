//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SPACESAVING_H
#define SPACESAVING_H

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

void SpaceSaving(char *p, off_t size, int k);

void WrapperSpaceSaving(char *p, off_t size, int kmer, double k, double l);

#endif // SPACESAVING_H

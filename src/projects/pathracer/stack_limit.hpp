
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <sys/time.h>
#include <sys/resource.h>
#include <errno.h>
#include <assert.h>

rlim_t stack_limit() {
   struct rlimit sl;
   int returnVal = getrlimit(RLIMIT_STACK, &sl);
   assert(returnVal != -1);
   return sl.rlim_cur;
}

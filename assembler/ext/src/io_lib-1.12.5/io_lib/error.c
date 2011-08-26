/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <sys/types.h>
#include <stdarg.h>

/*
 * Usage:
 *
 * errout(format, args...);
 */
void errout(char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);

    va_end(args);
}

/*
 * memmove() does not exist on SunOS 4.x, despite being an ANSI library call.
 *
 * void *memmove(void *to, const void *from, size_t len) {
 *    bcopy(from, to, len);
 *    return to;
 * }
 */

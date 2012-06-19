/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.h,v 1.9 2007/04/24 12:16:41 zerbino Exp $
 *
 */

#ifndef _FIB_H_
#define _FIB_H_

#include "globals.h"

typedef Coordinate(*voidcmp) (void *, void *);

/* functions for key heaps */
FibHeap *fh_makekeyheap(void);
FibHeapNode *fh_insertkey(FibHeap *, Coordinate, void *);
Coordinate fh_minkey(FibHeap *);
Coordinate fh_replacekey(FibHeap *, FibHeapNode *, Coordinate);
void *fh_replacekeydata(FibHeap *, FibHeapNode *, Coordinate, void *);

/* functions for void * heaps */
FibHeap *fh_makeheap(void);
voidcmp fh_setcmp(FibHeap *, voidcmp);
void *fh_setneginf(FibHeap *, void *);
FibHeapNode *fh_insert(FibHeap *, void *);

/* shared functions */
void *fh_extractmin(FibHeap *);
void *fh_min(FibHeap *);
void *fh_replacedata(FibHeapNode *, void *);
void *fh_delete(FibHeap *, FibHeapNode *);
void fh_deleteheap(FibHeap *);
FibHeap *fh_union(FibHeap *, FibHeap *);

#endif				/* _FIB_H_ */

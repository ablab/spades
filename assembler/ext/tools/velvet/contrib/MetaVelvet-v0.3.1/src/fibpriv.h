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
 * Copyright 1997, 1999-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without

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
 *	$Id: fibpriv.h,v 1.10 2007/10/09 09:56:46 zerbino Exp $
 *
 */

#ifndef _FIBPRIV_H_
#define _FIBPRIV_H_

#include "globals.h"

/*
 * specific node operations
 */
struct fibheap_el {
	FibHeapNode *fhe_p;
	FibHeapNode *fhe_child;
	FibHeapNode *fhe_left;
	FibHeapNode *fhe_right;
	void *fhe_data;
	Coordinate fhe_key;
	int fhe_degree;
	boolean fhe_mark;
};

static FibHeapNode *fhe_newelem(struct fibheap *);
static void fhe_initelem(FibHeapNode *);
static void fhe_insertafter(FibHeapNode * a, FibHeapNode * b);
static inline void fhe_insertbefore(FibHeapNode * a, FibHeapNode * b);
static FibHeapNode *fhe_remove(FibHeapNode * a);

/*
 * global heap operations
 */
struct fibheap {
	Coordinate(*fh_cmp_fnct) (void *, void *);
	RecycleBin *nodeMemory;
	IDnum fh_n;
	IDnum fh_Dl;
	FibHeapNode **fh_cons;
	FibHeapNode *fh_min;
	FibHeapNode *fh_root;
	void *fh_neginf;
	boolean fh_keys:1;
};

static void fh_initheap(FibHeap *);
static void fh_insertrootlist(FibHeap *, FibHeapNode *);
static void fh_removerootlist(FibHeap *, FibHeapNode *);
static void fh_consolidate(FibHeap *);
static void fh_heaplink(FibHeap * h, FibHeapNode * y, FibHeapNode * x);
static void fh_cut(FibHeap *, FibHeapNode *, FibHeapNode *);
static void fh_cascading_cut(FibHeap *, FibHeapNode *);
static FibHeapNode *fh_extractminel(FibHeap *);
static void fh_checkcons(FibHeap * h);
static void fh_destroyheap(FibHeap * h);
static int fh_compare(FibHeap * h, FibHeapNode * a, FibHeapNode * b);
static int fh_comparedata(FibHeap * h, Coordinate key,
			  void *data, FibHeapNode * b);
static void fh_insertel(FibHeap * h, FibHeapNode * x);

/*
 * general functions
 */
static inline IDnum ceillog2(IDnum a);

#endif				/* _FIBPRIV_H_ */

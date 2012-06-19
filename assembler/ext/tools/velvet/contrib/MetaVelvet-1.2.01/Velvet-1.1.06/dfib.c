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
 * Copyright 1997-2003 John-Mark Gurney.
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
 *	$Id: dfib.c,v 1.12 2007/10/19 13:09:26 zerbino Exp $
 *
 */
#include <limits.h>
#include <stdlib.h>

#include "recycleBin.h"
#include "dfib.h"

#include "dfibpriv.h"

#define BLOCKSIZE 10000
static DFibHeapNode *allocateDFibHeapNode(DFibHeap * heap)
{
  return (DFibHeapNode*)allocatePointer(heap->nodeMemory);
}

static void deallocateDFibHeapNode(DFibHeapNode * a, DFibHeap * heap)
{
	deallocatePointer(heap->nodeMemory, a);
}

IDnum dfibheap_getSize(DFibHeap * heap)
{
	return heap->dfh_n;
}

#define swap(type, a, b)		\
		do {			\
			type c;		\
			c = a;		\
			a = b;		\
			b = c;		\
		} while (0)		\

#define INT_BITS        (sizeof(IDnum) * 8)

static inline IDnum ceillog2(IDnum a)
{
	IDnum oa;
	IDnum i;
	IDnum b;
	IDnum cons;

	oa = a;
	b = INT_BITS / 2;
	i = 0;
	while (b) {
		i = (i << 1);
		cons = ((IDnum) 1) << b;
		if (a >= cons) {
			a /= cons;
			i = i | 1;
		} else
			a &= cons - 1;
		b /= 2;
	}
	if ((((IDnum) 1 << i)) == oa)
		return i;
	else
		return i + 1;
}

/*
 * Public Heap Functions
 */
DFibHeap *dfh_makekeyheap()
{
	DFibHeap *n;

	if ((n = (DFibHeap*)malloc(sizeof *n)) == NULL)
		return NULL;

	n->nodeMemory = newRecycleBin(sizeof(DFibHeapNode), BLOCKSIZE);
	n->dfh_n = 0;
	n->dfh_Dl = -1;
	n->dfh_cons = NULL;
	n->dfh_min = NULL;
	n->dfh_root = NULL;

	return n;
}

void dfh_deleteheap(DFibHeap * h)
{
	destroyRecycleBin(h->nodeMemory);
	if (h->dfh_cons != NULL)
		free(h->dfh_cons);
	free(h);
}

/*
 * Public Key Heap Functions
 */
DFibHeapNode *dfh_insertkey(DFibHeap * h, Time key, void *data)
{
	DFibHeapNode *x;

	if ((x = dfhe_newelem(h)) == NULL)
		return NULL;

	/* just insert on root list, and make sure it's not the new min */
	x->dfhe_data = data;
	x->dfhe_key = key;

	dfh_insertel(h, x);

	return x;
}

Time dfh_replacekey(DFibHeap * h, DFibHeapNode * x, Time key)
{
	Time ret;

	ret = x->dfhe_key;
	(void) dfh_replacekeydata(h, x, key, x->dfhe_data);

	return ret;
}

void *dfh_replacekeydata(DFibHeap * h, DFibHeapNode * x,
			 Time key, void *data)
{
	void *odata;
	Time okey;
	DFibHeapNode *y;
	int r;

	odata = x->dfhe_data;
	okey = x->dfhe_key;

	/*
	 * we can increase a key by deleting and reinserting, that
	 * requires O(lgn) time.
	 */
	if ((r = dfh_comparedata(h, key, data, x)) > 0) {
		/* XXX - bad code! */
		abort();
	}

	x->dfhe_data = data;
	x->dfhe_key = key;

	/* because they are equal, we don't have to do anything */
	if (r == 0)
		return odata;

	y = x->dfhe_p;

	if (okey == key)
		return odata;

	if (y != NULL && dfh_compare(h, x, y) <= 0) {
		dfh_cut(h, x, y);
		dfh_cascading_cut(h, y);
	}

	/*
	 * the = is so that the call from dfh_delete will delete the proper
	 * element.
	 */
	if (dfh_compare(h, x, h->dfh_min) <= 0)
		h->dfh_min = x;

	return odata;
}

/*
 * Public void * Heap Functions
 */
/*
 * this will return these values:
 *	NULL	failed for some reason
 *	ptr	token to use for manipulation of data
 */
void *dfh_extractmin(DFibHeap * h)
{
	DFibHeapNode *z;
	void *ret;

	ret = NULL;

	if (h->dfh_min != NULL) {
		z = dfh_extractminel(h);
		ret = z->dfhe_data;
		deallocateDFibHeapNode(z, h);
	}

	return ret;
}

void *dfh_replacedata(DFibHeapNode * x, void *data)
{
	void *odata = x->dfhe_data;
	x->dfhe_data = data;
	return odata;
}

void *dfh_delete(DFibHeap * h, DFibHeapNode * x)
{
	void *k;

	k = x->dfhe_data;
	dfh_replacekey(h, x, INT_MIN);
	dfh_extractmin(h);

	return k;
}

/*
 * begin of private element fuctions
 */
static DFibHeapNode *dfh_extractminel(DFibHeap * h)
{
	DFibHeapNode *ret;
	DFibHeapNode *x, *y, *orig;

	ret = h->dfh_min;

	orig = NULL;
	/* put all the children on the root list */
	/* for true consistancy, we should use dfhe_remove */
	for (x = ret->dfhe_child; x != orig && x != NULL;) {
		if (orig == NULL)
			orig = x;
		y = x->dfhe_right;
		x->dfhe_p = NULL;
		dfh_insertrootlist(h, x);
		x = y;
	}
	/* remove minimum from root list */
	dfh_removerootlist(h, ret);
	h->dfh_n--;

	/* if we aren't empty, consolidate the heap */
	if (h->dfh_n == 0)
		h->dfh_min = NULL;
	else {
		h->dfh_min = ret->dfhe_right;
		dfh_consolidate(h);
	}

	return ret;
}

static void dfh_insertrootlist(DFibHeap * h, DFibHeapNode * x)
{
	if (h->dfh_root == NULL) {
		h->dfh_root = x;
		x->dfhe_left = x;
		x->dfhe_right = x;
		return;
	}

	dfhe_insertafter(h->dfh_root, x);
}

static void dfh_removerootlist(DFibHeap * h, DFibHeapNode * x)
{
	if (x->dfhe_left == x)
		h->dfh_root = NULL;
	else
		h->dfh_root = dfhe_remove(x);
}

static void dfh_consolidate(DFibHeap * h)
{
	DFibHeapNode **a;
	DFibHeapNode *w;
	DFibHeapNode *y;
	DFibHeapNode *x;
	IDnum i;
	IDnum d;
	IDnum D;

	dfh_checkcons(h);

	/* assign a the value of h->dfh_cons so I don't have to rewrite code */
	D = h->dfh_Dl + 1;
	a = h->dfh_cons;

	for (i = 0; i < D; i++)
		a[i] = NULL;

	while ((w = h->dfh_root) != NULL) {
		x = w;
		dfh_removerootlist(h, w);
		d = x->dfhe_degree;
		/* XXX - assert that d < D */
		while (a[d] != NULL) {
			y = a[d];
			if (dfh_compare(h, x, y) > 0)
				swap(DFibHeapNode *, x, y);
			dfh_heaplink(h, y, x);
			a[d] = NULL;
			d++;
		}
		a[d] = x;
	}
	h->dfh_min = NULL;
	for (i = 0; i < D; i++)
		if (a[i] != NULL) {
			dfh_insertrootlist(h, a[i]);
			if (h->dfh_min == NULL
			    || dfh_compare(h, a[i], h->dfh_min) < 0)
				h->dfh_min = a[i];
		}
}

static void dfh_heaplink(DFibHeap * h, DFibHeapNode * y, DFibHeapNode * x)
{
	/* make y a child of x */
	if (x->dfhe_child == NULL)
		x->dfhe_child = y;
	else
		dfhe_insertbefore(x->dfhe_child, y);
	y->dfhe_p = x;
	x->dfhe_degree++;
	y->dfhe_mark = 0;
}

static void dfh_cut(DFibHeap * h, DFibHeapNode * x, DFibHeapNode * y)
{
	dfhe_remove(x);
	y->dfhe_degree--;
	dfh_insertrootlist(h, x);
	x->dfhe_p = NULL;
	x->dfhe_mark = 0;
}

static void dfh_cascading_cut(DFibHeap * h, DFibHeapNode * y)
{
	DFibHeapNode *z;

	while ((z = y->dfhe_p) != NULL) {
		if (y->dfhe_mark == 0) {
			y->dfhe_mark = 1;
			return;
		} else {
			dfh_cut(h, y, z);
			y = z;
		}
	}
}

/*
 * begining of handling elements of dfibheap
 */
static DFibHeapNode *dfhe_newelem(DFibHeap * h)
{
	DFibHeapNode *e;

	if ((e = allocateDFibHeapNode(h)) == NULL)
		return NULL;

	e->dfhe_degree = 0;
	e->dfhe_mark = 0;
	e->dfhe_p = NULL;
	e->dfhe_child = NULL;
	e->dfhe_left = e;
	e->dfhe_right = e;
	e->dfhe_data = NULL;

	return e;
}

static void dfhe_insertafter(DFibHeapNode * a, DFibHeapNode * b)
{
	if (a == a->dfhe_right) {
		a->dfhe_right = b;
		a->dfhe_left = b;
		b->dfhe_right = a;
		b->dfhe_left = a;
	} else {
		b->dfhe_right = a->dfhe_right;
		a->dfhe_right->dfhe_left = b;
		a->dfhe_right = b;
		b->dfhe_left = a;
	}
}

static inline void dfhe_insertbefore(DFibHeapNode * a, DFibHeapNode * b)
{
	dfhe_insertafter(a->dfhe_left, b);
}

static DFibHeapNode *dfhe_remove(DFibHeapNode * x)
{
	DFibHeapNode *ret;

	if (x == x->dfhe_left)
		ret = NULL;
	else
		ret = x->dfhe_left;

	/* fix the parent pointer */
	if (x->dfhe_p != NULL && x->dfhe_p->dfhe_child == x)
		x->dfhe_p->dfhe_child = ret;

	x->dfhe_right->dfhe_left = x->dfhe_left;
	x->dfhe_left->dfhe_right = x->dfhe_right;

	/* clear out hanging pointers */
	x->dfhe_p = NULL;
	x->dfhe_left = x;
	x->dfhe_right = x;

	return ret;
}

static void dfh_checkcons(DFibHeap * h)
{
	IDnum oDl;

	/* make sure we have enough memory allocated to "reorganize" */
	if (h->dfh_Dl == -1 || h->dfh_n > (1 << h->dfh_Dl)) {
		oDl = h->dfh_Dl;
		if ((h->dfh_Dl = ceillog2(h->dfh_n) + 1) < 8)
			h->dfh_Dl = 8;
		if (oDl != h->dfh_Dl)
			h->dfh_cons =
			    (DFibHeapNode **) realloc(h->dfh_cons,
						      sizeof *h->
						      dfh_cons *
						      (h->dfh_Dl + 1));
		if (h->dfh_cons == NULL)
			abort();
	}
}

static int dfh_compare(DFibHeap * h, DFibHeapNode * a, DFibHeapNode * b)
{
	if (a->dfhe_key < b->dfhe_key)
		return -1;
	if (a->dfhe_key == b->dfhe_key)
		return 0;
	return 1;
}

static int
dfh_comparedata(DFibHeap * h, Time key, void *data, DFibHeapNode * b)
{
	DFibHeapNode a;

	a.dfhe_key = key;
	a.dfhe_data = data;

	return dfh_compare(h, &a, b);
}

static void dfh_insertel(DFibHeap * h, DFibHeapNode * x)
{
	dfh_insertrootlist(h, x);

	if (h->dfh_min == NULL || x->dfhe_key < h->dfh_min->dfhe_key)
		h->dfh_min = x;

	h->dfh_n++;
}

Time dfibheap_el_getKey(DFibHeapNode * node)
{
	return node->dfhe_key;
}

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
 *	$Id: fib.c,v 1.10 2007/10/19 13:09:26 zerbino Exp $
 *
 */
#include <limits.h>
#include <stdlib.h>

#include "fib.h"
#include "recycleBin.h"

#include "fibpriv.h"

#define BLOCKSIZE 10000

static FibHeapNode *allocateFibHeapEl(FibHeap * heap)
{
	return allocatePointer(heap->nodeMemory);
}

static void deallocateFibHeapEl(FibHeapNode * a, FibHeap * heap)
{
	deallocatePointer(heap->nodeMemory, a);
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
 * Private Heap Functions
 */
static void fh_initheap(FibHeap * new)
{
	new->fh_cmp_fnct = NULL;
	new->nodeMemory = newRecycleBin(sizeof(FibHeapNode), BLOCKSIZE);
	new->fh_neginf = NULL;
	new->fh_n = 0;
	new->fh_Dl = -1;
	new->fh_cons = NULL;
	new->fh_min = NULL;
	new->fh_root = NULL;
	new->fh_keys = 0;
}

static void fh_destroyheap(FibHeap * h)
{
	h->fh_cmp_fnct = NULL;
	h->fh_neginf = NULL;
	if (h->fh_cons != NULL)
		free(h->fh_cons);
	h->fh_cons = NULL;
	free(h);
}

/*
 * Public Heap Functions
 */
FibHeap *fh_makekeyheap()
{
	FibHeap *n;

	if ((n = malloc(sizeof *n)) == NULL)
		return NULL;

	fh_initheap(n);
	n->fh_keys = 1;

	return n;
}

FibHeap *fh_makeheap()
{
	FibHeap *n;

	if ((n = malloc(sizeof *n)) == NULL)
		return NULL;

	fh_initheap(n);

	return n;
}

voidcmp fh_setcmp(FibHeap * h, voidcmp fnct)
{
	voidcmp oldfnct;

	oldfnct = h->fh_cmp_fnct;
	h->fh_cmp_fnct = fnct;

	return oldfnct;
}

void *fh_setneginf(FibHeap * h, void *data)
{
	void *old;

	old = h->fh_neginf;
	h->fh_neginf = data;

	return old;
}

FibHeap *fh_union(FibHeap * ha, FibHeap * hb)
{
	FibHeapNode *x;

	if (ha->fh_root == NULL || hb->fh_root == NULL) {
		/* either one or both are empty */
		if (ha->fh_root == NULL) {
			fh_destroyheap(ha);
			return hb;
		} else {
			fh_destroyheap(hb);
			return ha;
		}
	}
	ha->fh_root->fhe_left->fhe_right = hb->fh_root;
	hb->fh_root->fhe_left->fhe_right = ha->fh_root;
	x = ha->fh_root->fhe_left;
	ha->fh_root->fhe_left = hb->fh_root->fhe_left;
	hb->fh_root->fhe_left = x;
	ha->fh_n += hb->fh_n;
	/*
	 * we probably should also keep stats on number of unions
	 */

	/* set fh_min if necessary */
	if (fh_compare(ha, hb->fh_min, ha->fh_min) < 0)
		ha->fh_min = hb->fh_min;

	fh_destroyheap(hb);
	return ha;
}

void fh_deleteheap(FibHeap * h)
{
	destroyRecycleBin(h->nodeMemory);
	fh_destroyheap(h);
}

/*
 * Public Key Heap Functions
 */
FibHeapNode *fh_insertkey(FibHeap * h, Coordinate key, void *data)
{
	FibHeapNode *x;

	if ((x = fhe_newelem(h)) == NULL)
		return NULL;

	/* just insert on root list, and make sure it's not the new min */
	x->fhe_data = data;
	x->fhe_key = key;

	fh_insertel(h, x);

	return x;
}

Coordinate fh_minkey(FibHeap * h)
{
	if (h->fh_min == NULL)
		return (Coordinate) INT_MIN;
	return h->fh_min->fhe_key;
}

Coordinate fh_replacekey(FibHeap * h, FibHeapNode * x, Coordinate key)
{
	Coordinate ret;

	ret = x->fhe_key;
	(void) fh_replacekeydata(h, x, key, x->fhe_data);

	return ret;
}

void *fh_replacekeydata(FibHeap * h, FibHeapNode * x,
			Coordinate key, void *data)
{
	void *odata;
	Coordinate okey;
	FibHeapNode *y;
	int r;

	odata = x->fhe_data;
	okey = x->fhe_key;

	/*
	 * we can increase a key by deleting and reinserting, that
	 * requires O(lgn) time.
	 */
	if ((r = fh_comparedata(h, key, data, x)) > 0) {
		/* XXX - bad code! */
		abort();
	}

	x->fhe_data = data;
	x->fhe_key = key;

	/* because they are equal, we don't have to do anything */
	if (r == 0)
		return odata;

	y = x->fhe_p;

	if (h->fh_keys && okey == key)
		return odata;

	if (y != NULL && fh_compare(h, x, y) <= 0) {
		fh_cut(h, x, y);
		fh_cascading_cut(h, y);
	}

	/*
	 * the = is so that the call from fh_delete will delete the proper
	 * element.
	 */
	if (fh_compare(h, x, h->fh_min) <= 0)
		h->fh_min = x;

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
FibHeapNode *fh_insert(FibHeap * h, void *data)
{
	FibHeapNode *x;

	if ((x = fhe_newelem(h)) == NULL)
		return NULL;

	/* just insert on root list, and make sure it's not the new min */
	x->fhe_data = data;

	fh_insertel(h, x);

	return x;
}

void *fh_min(FibHeap * h)
{
	if (h->fh_min == NULL)
		return NULL;
	return h->fh_min->fhe_data;
}

void *fh_extractmin(FibHeap * h)
{
	FibHeapNode *z;
	void *ret;

	ret = NULL;

	if (h->fh_min != NULL) {
		z = fh_extractminel(h);
		ret = z->fhe_data;
#ifndef NO_FREE
		deallocateFibHeapEl(z, h);
#endif

	}

	return ret;
}

void *fh_replacedata(FibHeapNode * x, void *data)
{
	void *odata = x->fhe_data;
	x->fhe_data = data;
	return odata;
}

void *fh_delete(FibHeap * h, FibHeapNode * x)
{
	void *k;

	k = x->fhe_data;
	if (!h->fh_keys)
		fh_replacedata(x, h->fh_neginf);
	else
		fh_replacekey(h, x, (Coordinate) INT_MIN);
	fh_extractmin(h);

	return k;
}

/*
 * begin of private element fuctions
 */
static FibHeapNode *fh_extractminel(FibHeap * h)
{
	FibHeapNode *ret;
	FibHeapNode *x, *y, *orig;

	ret = h->fh_min;

	orig = NULL;
	/* put all the children on the root list */
	/* for true consistancy, we should use fhe_remove */
	for (x = ret->fhe_child; x != orig && x != NULL;) {
		if (orig == NULL)
			orig = x;
		y = x->fhe_right;
		x->fhe_p = NULL;
		fh_insertrootlist(h, x);
		x = y;
	}
	/* remove minimum from root list */
	fh_removerootlist(h, ret);
	h->fh_n--;

	/* if we aren't empty, consolidate the heap */
	if (h->fh_n == 0)
		h->fh_min = NULL;
	else {
		h->fh_min = ret->fhe_right;
		fh_consolidate(h);
	}

	return ret;
}

static void fh_insertrootlist(FibHeap * h, FibHeapNode * x)
{
	if (h->fh_root == NULL) {
		h->fh_root = x;
		x->fhe_left = x;
		x->fhe_right = x;
		return;
	}

	fhe_insertafter(h->fh_root, x);
}

static void fh_removerootlist(FibHeap * h, FibHeapNode * x)
{
	if (x->fhe_left == x)
		h->fh_root = NULL;
	else
		h->fh_root = fhe_remove(x);
}

static void fh_consolidate(FibHeap * h)
{
	FibHeapNode **a;
	FibHeapNode *w;
	FibHeapNode *y;
	FibHeapNode *x;
	IDnum i;
	IDnum d;
	IDnum D;

	fh_checkcons(h);

	/* assign a the value of h->fh_cons so I don't have to rewrite code */
	D = h->fh_Dl + 1;
	a = h->fh_cons;

	for (i = 0; i < D; i++)
		a[i] = NULL;

	while ((w = h->fh_root) != NULL) {
		x = w;
		fh_removerootlist(h, w);
		d = x->fhe_degree;
		/* XXX - assert that d < D */
		while (a[d] != NULL) {
			y = a[d];
			if (fh_compare(h, x, y) > 0)
				swap(FibHeapNode *, x, y);
			fh_heaplink(h, y, x);
			a[d] = NULL;
			d++;
		}
		a[d] = x;
	}
	h->fh_min = NULL;
	for (i = 0; i < D; i++)
		if (a[i] != NULL) {
			fh_insertrootlist(h, a[i]);
			if (h->fh_min == NULL
			    || fh_compare(h, a[i], h->fh_min) < 0)
				h->fh_min = a[i];
		}
}

static void fh_heaplink(FibHeap * h, FibHeapNode * y, FibHeapNode * x)
{
	/* make y a child of x */
	if (x->fhe_child == NULL)
		x->fhe_child = y;
	else
		fhe_insertbefore(x->fhe_child, y);
	y->fhe_p = x;
	x->fhe_degree++;
	y->fhe_mark = 0;
}

static void fh_cut(FibHeap * h, FibHeapNode * x, FibHeapNode * y)
{
	fhe_remove(x);
	y->fhe_degree--;
	fh_insertrootlist(h, x);
	x->fhe_p = NULL;
	x->fhe_mark = 0;
}

static void fh_cascading_cut(FibHeap * h, FibHeapNode * y)
{
	FibHeapNode *z;

	while ((z = y->fhe_p) != NULL) {
		if (y->fhe_mark == 0) {
			y->fhe_mark = 1;
			return;
		} else {
			fh_cut(h, y, z);
			y = z;
		}
	}
}

/*
 * begining of handling elements of fibheap
 */
static FibHeapNode *fhe_newelem(FibHeap * h)
{
	FibHeapNode *e;

	if ((e = allocateFibHeapEl(h)) == NULL)
		return NULL;

	fhe_initelem(e);

	return e;
}

static void fhe_initelem(FibHeapNode * e)
{
	e->fhe_degree = 0;
	e->fhe_mark = 0;
	e->fhe_p = NULL;
	e->fhe_child = NULL;
	e->fhe_left = e;
	e->fhe_right = e;
	e->fhe_data = NULL;
}

static void fhe_insertafter(FibHeapNode * a, FibHeapNode * b)
{
	if (a == a->fhe_right) {
		a->fhe_right = b;
		a->fhe_left = b;
		b->fhe_right = a;
		b->fhe_left = a;
	} else {
		b->fhe_right = a->fhe_right;
		a->fhe_right->fhe_left = b;
		a->fhe_right = b;
		b->fhe_left = a;
	}
}

static inline void fhe_insertbefore(FibHeapNode * a, FibHeapNode * b)
{
	fhe_insertafter(a->fhe_left, b);
}

static FibHeapNode *fhe_remove(FibHeapNode * x)
{
	FibHeapNode *ret;

	if (x == x->fhe_left)
		ret = NULL;
	else
		ret = x->fhe_left;

	/* fix the parent pointer */
	if (x->fhe_p != NULL && x->fhe_p->fhe_child == x)
		x->fhe_p->fhe_child = ret;

	x->fhe_right->fhe_left = x->fhe_left;
	x->fhe_left->fhe_right = x->fhe_right;

	/* clear out hanging pointers */
	x->fhe_p = NULL;
	x->fhe_left = x;
	x->fhe_right = x;

	return ret;
}

static void fh_checkcons(FibHeap * h)
{
	IDnum oDl;

	/* make sure we have enough memory allocated to "reorganize" */
	if (h->fh_Dl == -1 || h->fh_n > (1 << h->fh_Dl)) {
		oDl = h->fh_Dl;
		if ((h->fh_Dl = ceillog2(h->fh_n) + 1) < 8)
			h->fh_Dl = 8;
		if (oDl != h->fh_Dl)
			h->fh_cons =
			    (FibHeapNode **) realloc(h->fh_cons,
						     sizeof *h->
						     fh_cons *
						     (h->fh_Dl + 1));
		if (h->fh_cons == NULL)
			abort();
	}
}

static int fh_compare(FibHeap * h, FibHeapNode * a, FibHeapNode * b)
{
	if (a->fhe_key < b->fhe_key)
		return -1;
	if (a->fhe_key == b->fhe_key)
		return 0;
	return 1;
}

static int
fh_comparedata(FibHeap * h, Coordinate key, void *data, FibHeapNode * b)
{
	FibHeapNode a;

	a.fhe_key = key;
	a.fhe_data = data;

	return fh_compare(h, &a, b);
}

static void fh_insertel(FibHeap * h, FibHeapNode * x)
{
	fh_insertrootlist(h, x);

	if (h->fh_min == NULL
	    || (h->fh_keys ? x->fhe_key <
		h->fh_min->fhe_key : h->fh_cmp_fnct(x->fhe_data,
						    h->fh_min->fhe_data) <
		0))
		h->fh_min = x;

	h->fh_n++;
}

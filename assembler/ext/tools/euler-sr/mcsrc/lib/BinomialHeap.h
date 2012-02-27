/***************************************************************************
 * Title:          BinomialHeap.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BINOMIAL_HEAP_H_
#define BINOMIAL_HEAP_H_
#include <assert.h>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

template<typename K, typename V>
class BinomialHeapNode {
 public:
	K key;
	V value;
	BinomialHeapNode<K,V> *parent, *left, *sibling, *preSibling, **extRef;
	ssize_t degree;
	BinomialHeapNode(K k, V &v) {
		parent = left = sibling = preSibling = NULL;
		extRef = NULL;
		key = k;
		value = v;
		degree = 0;
		
	}
	BinomialHeapNode() {
		parent = left = sibling = NULL;
	}
};

template<typename K, typename V>
class BinomialHeap {
 public:
	ssize_t nNodes;
	BinomialHeapNode<K,V> *head;
	BinomialHeap() {
		head = NULL;
		nNodes = 0;
	}
	BinomialHeapNode<K,V> *MakeBinomialHeapNode(K key, V &value) {
		BinomialHeapNode<K,V> *node = new BinomialHeapNode<K,V>(key, value);
		return node;
	}

	BinomialHeapNode<K,V> *MakeBinomialHeap() {
		BinomialHeapNode<K,V> H = new BinomialHeapNode<K,V>;
		H->head = NULL; // should already be done in the constructor.
		return H;
	}

	ssize_t Empty() {
		return head == NULL;
	}

	BinomialHeapNode<K,V> HeapMinimum() {
		BinomialHeapNode<K,V> *x, *y;
		y = NULL;
		x = head;
		K min;

		ssize_t first = 1;
		while (x != NULL) {
			if (first or x->key < min) {
				first = 0;
				min = x->key;
				y = x;
			}
			x = x->sibling;
		}
		return y;
	}

	void Link(BinomialHeapNode<K,V> *y,
						BinomialHeapNode<K,V> *z) {
		y->parent = z;
		y->sibling = z->left;
		z->left = y;
		z->degree++;
	}

	BinomialHeapNode<K, V> *Merge(BinomialHeapNode<K,V> *h1,
			             	BinomialHeapNode<K,V> *h2) {

		/*
		 * Merge the ordered lists of h1 and h2 into a single ordered list.
		 */

		BinomialHeapNode<K,V> *cur1, *cur2, *h, *cur, *min;

		cur1 = h1; cur2 = h2;
		h = NULL;
		cur = NULL;
		while (cur1 != NULL and 
			cur2 != NULL) {
			if (cur1->degree < cur2->degree) {
				min = cur1;
				cur1 = cur1->sibling;
			}
			else {
				min  = cur2;
				cur2 = cur2->sibling;
			}
			if (h == NULL) {
				// Create the head if this is the first time.
				h   = min;
				cur = h;
			}
			else {
				cur->sibling = min;
				cur = cur->sibling;
			}
		}
		// Append the rest of cur1 if it exists
		if (cur1 != NULL) {
			if (h == NULL) {
				cur = h = cur1;
			}
			else {
				cur->sibling = cur1;
			}
		}
		// 
		// Append the rest of cur2 if it exists
		if (cur2 != NULL) {
			if (h == NULL) {
				cur = h = cur2;
			}
			else {
				// append this element.
				cur->sibling = cur2;
			}
		}
		return h;
	}

	BinomialHeapNode<K,V> *Union(BinomialHeapNode<K,V> *h1,
															 BinomialHeapNode<K,V> *h2) {
		BinomialHeapNode<K,V>* h;
		h = Merge(h1, h2);
		if (h == NULL) {
			return h;
		}
		BinomialHeapNode<K,V> *x, *prevX = NULL, *nextX = NULL;
		x = h;
		assert(x != NULL);
		nextX = x->sibling;
		while (nextX != NULL) {
			if (x->degree != nextX->degree or
					(nextX->sibling != NULL and nextX->sibling->degree == x->degree)) {
				prevX = x;
				x = nextX;
			}
			else {
				if (x->key <= nextX->key) {
					x->sibling = nextX->sibling;
					Link(nextX, x);
				}
				else {
					if (prevX == NULL) {
						h = nextX;
					}
					else {
						prevX->sibling  = nextX;
					}
					Link(x, nextX);
					x = nextX;
				}
			}
			assert(x != NULL);
			nextX = x->sibling;
		}
		return h;
	}

	BinomialHeapNode<K,V>* Insert(K key, V value) {
		BinomialHeapNode<K,V> *hPrime;
		hPrime = MakeBinomialHeapNode(key,value);
		head = Union(head, hPrime);
		assert(head->parent == NULL);
		return hPrime;
	}

	//
	// Pull the lowest value node off the heap.
	// Return that node so we know what it is.
	//
	BinomialHeapNode<K, V> *ExtractMin() {

		BinomialHeapNode<K, V> *x, *min, *prevMin, *prevX;
		// Leave if the list is empty.
		if (head == NULL) {
			return NULL;
		}
		
		//
		// Find the minimum.
		//
		min = x = head;
		prevMin = prevX = NULL;
		while (x != NULL) {
			if (x->key <= min->key) {
				min = x;
				prevMin = prevX;
			}
			prevX = x;
			x = x->sibling;
		}
		//
		// Remove the minimum.
		//
		if (min == head) {
			// If the min is the head, must update that
			// instead of the prev.
			head = min->sibling;
		}
		else {
			// Min is a non-head node, make 
			// the prev skip past it.
			prevMin->sibling = min->sibling;
		}
		
		// 
		// Reverse the order of the children.
		//UNUSED// BinomialHeapNode<K, V> *reverse, *sibling;
		BinomialHeapNode<K, V> *cur, *next, *prev;
		cur  = min->left;
		if (cur != NULL) {
			next = cur->sibling;
			prev = NULL;
			while (cur->sibling != NULL) {
				
				// 
				// Record where this will go next
				next  = cur->sibling;

				// 
				// Reverse the current.
				cur->sibling = prev;
				cur->parent = NULL;
				// 
				// Record where the next reversal will be
				prev = cur;

				//
				// Advance
				cur  = next;
			}
			// handle the last case.
			cur->sibling = prev;
			cur->parent = NULL;
		}
		// Add the reversed list to the head.
		if (head != NULL)
			assert(head->parent == NULL);
		head = Union(head, cur);
		if (head != NULL)
			assert(head->parent == NULL);
		return min;
	}		

	void Exchange(BinomialHeapNode<K,V> *x, BinomialHeapNode<K,V> *y) {

		K tmpKey;
		V tmpValue;

		tmpKey   = x->key;
		tmpValue = x->value;

		x->key  = y->key;	x->value = y->value;
		y->key = tmpKey; y->value = tmpValue;
		// external reference is a reference to 
		if (x->extRef != NULL && y->extRef != NULL) {
			*(x->extRef) = y;
			*(y->extRef) = x;
			BinomialHeapNode<K,V> **xExtRefCopy;
			xExtRefCopy = x->extRef;
			x->extRef = y->extRef;
			y->extRef = xExtRefCopy;
		}
	}
	
	void Swap(BinomialHeapNode<K,V> *x, BinomialHeapNode<K,V> *y) {
		
		//
		// X is the parent of Y
		//
		assert(y->parent == x);
		
		BinomialHeapNode<K,V> *xParent, *yLeft, *preX, *preY, 
			*xSibling, *ySibling;
		
		preX = FindPreSibling(x);
		preY = FindPreSibling(y);
		xSibling = x->sibling;
		ySibling = y->sibling;

		//
		// Move the parent down.
		
		x->left = y->left;
		if (x->left != NULL) {
			AssignParent(x->left, x);
			x->parent = y;
		}
		
		//
		// Move the child up.
		y->parent = xParent;
		y->left   = x;
		if (xParent != NULL) {
			if (xParent->left == x) {
				xParent->left = y;
			}
		}

		//
		// Move y into x's old list.
		y->sibling = xSibling;
		if (preX != NULL)
			preX->sibling = y;

		//
		// Move x into y's old list.
		x->sibling = ySibling;
		if (preY != NULL)
			preY->sibling = x;
	}

	void AssignParent(BinomialHeapNode<K,V>* node, 
										BinomialHeapNode<K,V> *parent) {
		while (node != NULL) {
			node->parent = parent;
			node = node->sibling;
		}
	}

	BinomialHeapNode<K,V> *FindPreSibling(BinomialHeapNode<K,V> *x) {
		BinomialHeapNode<K,V> *parent, *preSibling;
		parent = x->parent;
		if (parent == NULL) {
			preSibling = head;
		}
		else {
			preSibling = parent->left;
		}

		while (preSibling != NULL) {
			if (preSibling->sibling == x)
				break;
			preSibling = preSibling->sibling;
		}
		return preSibling;
	}

	void DecreaseKey(BinomialHeapNode<K,V> *x, K key) {
		assert(key < x->key);
		x->key = key;
		BinomialHeapNode<K,V> *y, *z;
		y = x;
		z = y->parent;
		while (z != NULL and y->key < z->key) {

			//			BinomialHeapNode<K,V> tmp;
			Exchange(y,z);
			y = z;
			z = y->parent;
		}
	}

	void Delete(BinomialHeapNode<K,V> *x) {
		DecreaseKey(x, -9999999);
		ExtractMin(head);
	}
	
	void Print() {
		stringstream curLine;
		Print(curLine, head);
	}


	void Print2() {
		Print2(head, "");
	}

	void Print2(BinomialHeapNode<K,V> *start, string padding) {
		cout << "( ";
		BinomialHeapNode<K,V> *cur = start;
		if (start == NULL) {
			cout << endl << endl;
			return;
		}
		while (cur != NULL) {
			cout << padding << cur->key << " ";
			if (cur->left != NULL) {
				Print2(cur->left, padding);
			}
			cur = cur->sibling;
		}
		cout << ")";
		if (start == head)
			cout << endl;
		
	}
		
		

	void Print(stringstream &curLine, BinomialHeapNode<K,V> *cur) { 
		if (head == NULL) {
			cout << endl << endl;
			return;
		}

		BinomialHeapNode<K,V> *start = cur;
		
		// print this node.
		ssize_t prefixLength = curLine.str().size();
		curLine << "(" << cur->value << ", " << cur->key << ") ";
		// recursively print all of its children.
		if (cur->left != NULL) {
			Print(curLine, cur->left);
		}
		else {
			cout << curLine.str() << endl;
		}

		// Print all siblings below
		while (cur->sibling != NULL) {
			ssize_t i;
			stringstream siblingLine;
			for (i = 0; i < prefixLength ; i++) 
				siblingLine << " ";
			//			Print(siblingLine, cur->sibling);
			assert(cur != cur->sibling);
			cur = cur->sibling;
		}
		if (start == head) {
			cout << std::endl << std::endl;
		}
	}
		
	ssize_t Free() {
		return Free(head);
	}

	ssize_t Free(BinomialHeapNode<K,V> *cur) {
		ssize_t numFreed = 0;
		BinomialHeapNode<K,V> *next;
		while (cur != NULL) {
			numFreed += Free(cur->left);
			next = cur->sibling;
			++numFreed;
			delete cur;
			cur = next;
		}
		return numFreed;
	}
};
	

#endif

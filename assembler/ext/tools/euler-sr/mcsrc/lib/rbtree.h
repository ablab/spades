/***************************************************************************
 * Title:          rbtree.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _RBTREE_
#define _RBTREE_

#include <ostream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define red 0
#define black 1

template <typename T> inline
ssize_t Greater(T a, T b) {
  return a > b;
}

template <typename T> inline 
ssize_t Greater(T* a, T *b) {
  return *a > *b;
}


template <typename T> inline
ssize_t Less(T a, T b) {
  std::cout << "using rbtree less " << std::endl;
  return a < b;
}
template <typename T> inline
ssize_t Less(T *a, T *b) {
  return *a < *b;
}
template <typename T> inline
ssize_t Equals(T a, T b) {
  return a == b;
}

template <typename T> inline
ssize_t Equals(T *a, T *b) {
  return *a == *b;
}


template <typename T> 
class RBTreeNode {

public:
  T data;
  RBTreeNode<T> *leftChild;
  RBTreeNode<T> *rightChild;
  RBTreeNode<T> *parent;
  ssize_t color;
  RBTreeNode() {
    leftChild = NULL;
    rightChild= NULL;
    parent    = NULL;
    color = red;
  }
  bool operator<(const T &other) const {
    return (this->data < other);
  }
  bool operator>(const T& other) const {
    return this->data > other;
  }
  bool operator==(const T& other) const {
    return this->data == other;
  }

};


template <typename T>
class RBTree {
  // red-black tree implementation, based off of method presented in CLR
  ssize_t _size;
public:
  RBTreeNode<T> *root;
  RBTreeNode<T> *nil; //sentinal used by rbdelete
  ssize_t size() { return _size; }
  RBTree() { 
    root = NULL;
    nil = new RBTreeNode<T>;
    nil->color = -1;
    _size = 0;
  }
  void RightRotate(RBTreeNode<T> *parent) {
    RBTreeNode<T>  *child;
    child = parent->leftChild;
    assert(child != NULL);
    // Perform the swap
    parent->leftChild = child->rightChild;
    assert(parent->leftChild != NULL);
    child->rightChild = parent;
    assert(child->rightChild != NULL);

    // Fix the parents of the swapped nodes
    child->parent  = parent->parent;
    parent->parent = child;

    // Fix the parents of the children of the swapped nodes.
    if (parent->leftChild != nil)
      parent->leftChild->parent = parent;
    if (child->rightChild != nil)
      child->rightChild->parent = child;


    if (child->parent != NULL && child->parent->leftChild == parent)
      child->parent->leftChild = child;
    else if (child->parent != NULL && child->parent->rightChild == parent)
      child->parent->rightChild = child;
    
    // have correct entry into the tree structure
    if (IsRoot(child))
      root = child;
  }

  void LeftRotate(RBTreeNode<T> *parent) {
    RBTreeNode<T> *child; 
    child = parent->rightChild;

    // Perform the swap
    parent->rightChild = child->leftChild;
    child->leftChild   = parent;

    // Fix the parents of the swapped nodes
    child->parent  = parent->parent;
    parent->parent = child;

    // Fix the parents of the subtrees.
    assert(parent->rightChild != NULL);
    assert(parent->leftChild != NULL);
    if (parent->rightChild != nil) 
      parent->rightChild->parent = parent;
    if (child->leftChild != NULL)
      child->leftChild->parent = child;

    // Fix entry point into this subtree
    // child->parent is now grandparent. 
    if (child->parent != NULL && child->parent->leftChild == parent)
      child->parent->leftChild = child;
    else if (child->parent != NULL && child->parent->rightChild == parent)
      child->parent->rightChild = child;
    
    if (IsRoot(child))
      root = child;
  }

  ssize_t IsRoot(RBTreeNode<T> *node) {
    return (node->parent == NULL);
  }

  ssize_t Color(RBTreeNode<T> *node) {
    if (node == NULL)
      return black;
    else 
      return node->color;
  }

  RBTreeNode<T>* Parent(RBTreeNode<T> *node) {
    if (node == NULL)
      return NULL;
    else
      return node->parent;
  }
  
  RBTreeNode<T>* Left(RBTreeNode<T> *node) {
    if (node == NULL) 
      return NULL;
    else
      return node->leftChild;
  }
  
  RBTreeNode<T>* Right(RBTreeNode<T> *node) {
    if (node == NULL)
      return NULL;
    else
      return node->rightChild;
  }

  RBTreeNode<T>* UnbalancedInsert(const T &data, 
				  RBTreeNode<T> *&node, 
				  RBTreeNode<T> *parent) {
    if (node == NULL || node == nil) {
      node = new RBTreeNode<T>;
      node->data = data;
      node->parent = parent;
      node->color = red;
      node->leftChild = nil;
      node->rightChild = nil;
      return node;
    }
    else {
			// TODO: Fix compiler warning "warning: '<anonymous>' may be used uninitialized in this function"
      if (Greater(data, node->data)) {
				return UnbalancedInsert(data, node->rightChild, node);
      } else {
				return UnbalancedInsert(data, node->leftChild, node);
			}
    }
  }
  
  RBTreeNode<T>* Insert(const T &data) { 
    _size++;
    RBTreeNode<T> *node;
    node = UnbalancedInsert(data, root, NULL);
    RBBalance(node);
    root->color = black;
#ifdef DEBUG
    CheckColoring(root);
#endif

    return node;
  }

  RBTreeNode<T>* RBBalance(RBTreeNode<T> *node) {
    RBTreeNode<T> *ynode;
		//UNUSED// RBTreeNode<T> *temp;
    while (!IsRoot(node) && Color(Parent(node)) == red) {
      if (Parent(node) == Left(Parent(Parent(node)))) {
	ynode = Right(Parent(Parent(node)));
	if (Color(ynode) == red) {
	  Parent(node)->color = black;
	  ynode->color = black;
	  Parent(Parent(node))->color = red;
	  node = Parent(Parent(node));
	}
	else {
	  if (node == Right(Parent(node))){ 
	    node = Parent(node);
	    LeftRotate(node);
	  }
	  Parent(node)->color = black;
	  Parent(Parent(node))->color = red;
	  RightRotate(Parent(Parent(node)));
	}
      }
      else {
	// Symmetric conditions
	ynode = Left(Parent(Parent(node)));
	if (Color(ynode) == red) {
	  Parent(node)->color = black;
	  ynode->color = black;
	  Parent(Parent(node))->color = red;
	  node = Parent(Parent(node));
	}
	else {
	  if (node == Left(Parent(node))){ 
	    node = Parent(node);
	    RightRotate(node);
	  }
	  Parent(node)->color = black;
	  Parent(Parent(node))->color = red;
	  LeftRotate(Parent(Parent(node)));
	}
      }      
    }
    return node;
  }

  void Print(std::ostream &out) {
    std::string padding = "";
    Print(out, root, padding);
  }

  void Print(std::ostream &out, RBTreeNode<T> *node, std::string padding = "") {

    if (node == NULL)
      node = root;

    if (root == NULL) return; // don't even try to print if there is nothing here

    out << padding << " " << node->color << " " << node->data << std::endl;
    
    if (node->leftChild != nil)   // only recurse on non-null nodes
      Print(out, node->leftChild, padding+" l ");
    else
      out << padding + " l " << "nil" << std::endl;
    if (node->rightChild != nil)
      Print(out, node->rightChild, padding + " r ");
    else
      out << padding + " r " << "nil" << std::endl;
  }


  void GetMaxDepth(ssize_t &maxDepth, RBTreeNode<T> *node= NULL , ssize_t curDepth=0) {
    if (node == NULL)
      node = root;
    if (node == NULL)
      return;

    if (curDepth > maxDepth)
      maxDepth = curDepth;
    
    if (node->leftChild != nil)
      GetMaxDepth(maxDepth, node->leftChild, curDepth++);
    
    if (node->rightChild != nil)
      GetMaxDepth(maxDepth, node->rightChild, curDepth++);
  }

  RBTreeNode<T>* Delete(RBTreeNode<T> *z) {
    ssize_t doFixup = 0;
    RBTreeNode<T> *y, *x;
    _size--;
    ssize_t wasRoot = 0;
    if (IsRoot(z)) {
      //      std::cout << "removing root reference " << z << std::endl;
      wasRoot = 1;
    }

    // find what node to detach from the graph.
    if (Left(z) == nil || Right(z) == nil) 
      // z has one child, detach z itself
      y = z;
    else
      // z has two children, detach the successor, but copy 
      // all information from the successor into z
      y = Successor(z);

    // Detach y from the graph.
    //    std::cout << "succ (" << z << " " << y << ")" << std::endl;
    if (Left(y) != nil) 
      x = Left(y);
    else
      x = Right(y);
    //    std::cout << "x: " << x << std::endl;
    x->parent = y->parent;

    if (IsRoot(y)) {
      root = x;
      // If nothing is left of the tree, record that as NULL, rather than nil
      if (root == nil)
	root = NULL;
    }
    else {
      // do the detaching, link past y
      if (y == Left(Parent(y)))
	y->parent->leftChild = x;
      else
	y->parent->rightChild = x;
    }
    
    if (y != z) {
      // Differ from CLR here.  This is the case that z had two children and 
      // as replaced by its successor.  CLR splices out the successor, and 
      // replaces z's data with y.  That's a problem when we have external 
      // references to the tree structure, so when z is deleted, it needs to go!
      // Rather than simply copying y into z, fully replace z with y. 
      // The problem is that we check to balance the tree if y's color is black,
      // not if z's color is black.  But y has become the new z here, so
      // we need to store a variable that says we need to do the fixup.

      if (y->color == black)
	doFixup = 1;
      else 
	doFixup = 0;

      // copy z's structure into y.  leave y's data alone
      y->leftChild  = z->leftChild;
      y->rightChild = z->rightChild;
      //      std::cout << "z's children: " << z->leftChild << " " << z->rightChild << std::endl;
      if (! Nil(y->leftChild) ) 
	y->leftChild->parent = y;

      if ( ! Nil(y->rightChild) ) 
	y->rightChild->parent = y;

      y->parent = z->parent;
      y->color  = z->color;
      if (IsRoot(z))
	root = y;
      else {
	if (Left(z->parent) == z)
	  z->parent->leftChild = y;
	if (Right(z->parent) == z)
	  z->parent->rightChild = y;
      }

      // fix problems with the pointer of x
      if (x->parent == z) 
	// z has been replaced by y, so make the parent of x y now
	x->parent = y;
    }
    else {
      if (y->color == black) 
	doFixup = 1;
    }
      
    /*    
	  std::cout << "rbtree before fixup: " << std::endl;
    */
    if (doFixup)
      RBDeleteFixup(x);

    /* 
       if (wasRoot) {
       std::cout << "old root: " << z << " new root " << root << std::endl;
       }
    */

    #ifdef DEBUG
    CheckColoring(root);
    #endif
    return y;
  }

  void RBDeleteFixup(RBTreeNode<T> *x) {
    RBTreeNode<T> *w;
    while (!IsRoot(x) && x->color== black) {
      if (x == Left(Parent(x))) { 
	w = Right(Parent(x));
	if (w->color == red) {
	  w->color = black;
	  Parent(x)->color = red;
	  LeftRotate(Parent(x));
	  w = Right(Parent(x));
	}
	if (Color(Left(w)) == black && Color(Right(w)) == black) {
	  w->color = red;
	  x = Parent(x);
	}
	else {
	  if (Right(w)->color == black) {
	    Left(w)->color = black;
	    w->color = red;
	    RightRotate(w);
	    w = Right(Parent(x));
	  }
	  w->color = Parent(x)->color;
	  Parent(x)->color = black;
	  Right(w)->color = black;
	  LeftRotate(Parent(x));
	  x = Root(x);
	}
      } // end x is left child
      else {
	w = Left(Parent(x));
	if (w->color == red) {
	  w->color = black;
	  Parent(x)->color = red;
	  RightRotate(Parent(x));
	  w = Left(Parent(x));
	}
	if (Color(Right(w)) == black && Color(Left(w)) == black) {
	  w->color = red;
	  x = Parent(x);
	}
	else {
	  if (Left(w)->color == black) {
	    Right(w)->color = black;
	    w->color = red;
	    LeftRotate(w);
            w = Left(Parent(x));
	  }
	  w->color = Parent(x)->color;
	  Parent(x)->color = black;
	  Left(w)->color = black;
	  RightRotate(Parent(x));
	  x = Root(x);
	}	
      }
    }
    x->color = black; // fix the color of the root
  }

  RBTreeNode<T> *Root(RBTreeNode<T> *node) {
    assert(node != NULL);
    while (node->parent != NULL)
      node = node->parent;
    return node;
  }

  ssize_t Nil(RBTreeNode<T> *node) {
    return (node == NULL || node == nil);
  }

  ssize_t Pop(T& data) {
    RBTreeNode<T>* min;
    min = TreeMinimum();
    
    if (min != NULL) {
      data = min->data;
      Delete(min);
      return 1;
    }
    // otherwise no min found
    return 0;
  }

  ssize_t PopMax(T& data) {
    RBTreeNode<T>* min;
    min = TreeMaximum();
    
    if (min != NULL) {
      data = min->data;
      Delete(min);
      return 1;
    }
    // otherwise no min found
    return 0;
  }

  RBTreeNode<T> *TreeMaximum(RBTreeNode<T> *node=NULL) {
    if (node == NULL)
      node = root;
    if (node == nil || node == NULL)
      return NULL;
    while (node->rightChild != nil) 
      node = node->rightChild;
    return node;
  }

  RBTreeNode<T> *TreeMinimum(RBTreeNode<T> *node=NULL) {
    if (node == NULL)
      node = root;
    if (node == nil || node == NULL)
      return NULL;
    while (node->leftChild != nil) 
      node = node->leftChild;
    return node;
  }

  RBTreeNode<T> *Successor(RBTreeNode<T>*node) {
    if (node->rightChild != nil)
      return TreeMinimum(node->rightChild);
    
    RBTreeNode<T> *tempNode;
    tempNode = node->parent;
    while (tempNode != nil && node == Right(tempNode)) {
      node = tempNode;
      tempNode = tempNode->parent;
    }
    return tempNode;
  }
	
	ssize_t FindData(const T &data) {
		if (Find(data) == NULL)
			return 0;
		else
			return 1;
	}
 
  RBTreeNode<T> *Find(const T &data, RBTreeNode<T>*node=NULL) {
    // do binary search for value;
    if (node == NULL) {
      node = root;
    }
    if (node == NULL || node == nil) 
      return NULL;

    if (Equals(node->data, data)) 
      return node;
    
    if (Greater(node->data, data))
      return Find(data, node->leftChild);
    
    if (Less(node->data, data))
      return Find(data, node->rightChild);
    
    assert(printf("Find made it through too many cases") == 0);
    return NULL;
  }

  ssize_t CheckColoring(RBTreeNode<T> *node) {
    ssize_t leftBlackCount, rightBlackCount;
    if (node == NULL)
      return 1;
    /*
      if (node == root) {
      std::cout << *this;
      }
    */
    if (node == nil) {
      return 1;
    }
    else {
      leftBlackCount = CheckColoring(node->leftChild);
      rightBlackCount = CheckColoring(node->rightChild);
      // If all is legit with the tree, left and righ
      // have the same number of black children
      if (leftBlackCount != rightBlackCount) {
	std::cout << *this;
      }
      assert(leftBlackCount == rightBlackCount);

      if (node->color == black) {
	return leftBlackCount + 1;
      }
      else {
	return leftBlackCount;
      }
    }
  }
  // This defines a find that should be subclassed by finds that are specific
  // to the application in mind.  Specifically, there will be a find 
  // that looks for the member edges that are contained inside any edge 
  // stored in the tree.
};

template<typename T>
std::ostream& operator<<(std::ostream &out, const RBTree<T> &tree) {
  tree.Print(out, tree.root);
  return out;
}
#endif

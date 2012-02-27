/***************************************************************************
 * Title:          BWdictionary.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/*
** Copyright 2004 Ross Lippert
** 
** This file is part of bbbwt.
** 
** bbbwt is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** bbbwt is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with bbbwt; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
#ifndef BWDICTIONARY_H
#define BWDICTIONARY_H

#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include "word.h"
#include "bitio.h"

namespace BW {

  // the vnode holds all of the query and modification methods
  // a dictionary is a special kind of vnode.
  // outside programs should create dictionaries, not vnodes
  template <class BitDict>
  class vnode_T {
    // red/black tree implementation
    // Memory is divided into bitdictionaries and nodes
    // bitdictionaries implement an interface similar to the dictionary, but
    //  are of fixed size.  They form the leaf set.
    // nodes implement the tree structure on top of the bitdictionaries.  Nodes
    //  store a minimum amount of data, two child pointers and two counts, the
    //  total bit and sum of bit of the left child (the right counts can be infered
    //  if the total count and sum for the node is known).
    // A vnode is something like an iterator/smart pointer.  It is a virtual node
    //  in the sense that it has data which is either a node or a bitdictionary
    //  and has extra data in it like total counts.
    // dictionary is just a vnode with a few extras.

    struct node;
    typedef union {
      node *u_node; 
      typename BitDict::bits *u_leaf;} node_or_leaf;
    
    struct node {
      // horrible bit tricks
			// HAS BITFIELD
      word_t rb:1;
      word_t sizeL:63;   // left size
      
      word_t L_is_leaf:1;
      word_t R_is_leaf:1;
      word_t sumL:62;    // left sum
    
      node_or_leaf L,R;     // lesser/greater pointers

      // nodes are created from leaves, and this is how
      node(typename BitDict::bits *l);
      // default node creation should hardly ever be used (is used in fscan)
      node() {}
    };

    // vnode members
  protected:
		// HAS BITFIELD
    word_t  _is_leaf:1;  // is my data bitdictionary
    word_t   _sum:63;     // sum of my data
    word_t  _size;       // number of bits of my data
    node_or_leaf data;
    vnode_T<BitDict> *parent;   // used to make a link list on the stack when needed

    // only dictionary's and vnodes will create and destroy vnodes
    ~vnode_T() {};
    vnode_T(): _is_leaf(1),_sum(0),_size(0),parent(0) {};
    
  public:
    inline bool is_leaf() const {return _is_leaf;}
    inline word_t size() const {return _size;}
    inline word_t  sum() const {return _sum;}

    // compute the bitvalue at position i
    void lookup(word_t pos,bool &b) const;
    // compute sum before i to sum
    void lookup(word_t pos,word_t &sum) const;
    // compute bitvalue at i and sum before i to sum
    void lookup(word_t pos,word_t &sum,bool &b) const;

    word_t find0(word_t) const; // find nth 0
    word_t find1(word_t) const; // find nth 1

    // insert bit value at position p.  also compute sum before p
    void insert(word_t p,bool b,word_t &s); // sum and insert bit value

    // IO routines to write to streams
    bool fprint(FILE *) const;
    bool fscan(FILE *) ;

    // an ascii dumper for debugging
    void fdump(FILE *) const;

    // reclaim lost memory in the bitdictionaries by merging them
    // do not call append after freeze unless thaw-ed.
    void freeze();

    // restore/check rb property to ensure it will be ok to do future
    // appends.
    void thaw();

    // useful for checking out the tree, and making sure it is ok
    bool   rb_check(word_t &) const;
    word_t  depth() const;
    word_t num_nodes() const;

  protected:
    // turns a vnode into the vnode which is logically to the left/right of it.
    //  assumes that the vnode points to a node and not a bitdictionary.
    //  adjusts all count data in the vnode appropriately.
    inline vnode_T<BitDict> &goleft() {
      _size = data.u_node->sizeL;
      _sum  = data.u_node->sumL;
      _is_leaf  = data.u_node->L_is_leaf;
      data  = data.u_node->L;
      return (*this);
    }
    inline vnode_T<BitDict> &goright() {
      _size -= data.u_node->sizeL;
      _sum  -= data.u_node->sumL;
      _is_leaf  = data.u_node->R_is_leaf;
      data  = data.u_node->R;
      return (*this);
    }

    // delete all pointers (nodes and bitdictionaries) under this vnode
    void destroy();

  private:
    // internal print function for outputting bitstreams
    bool fprint(bitio &) const;
    // internal scan to read from bitstreams
    bool fscan(bitio &);

    // turn the tree into that tree where every L pointer points to
    //  a bit dictionary, and every R but one points to a node
    //  and the depth first order of the bitdictionaries is preserved.
    //  Counts are adjusted.
    //  Basically turns the tree into a link list.
    node *flatten(word_t &szR,word_t &smR);

    // a very poor way to impose rb structure on a linked list in R
    bool regrow(bool rb);

    // only run on a flattened dictionary!  Causes the bitdictionaries
    //  which are adjacent in the linked list to be merged, with resulting
    //  empty bitdictionaries/nodes automatically deleted.
    void merge();

    // a very good way to impose rb structure on a linked list in R.
    //  rebuilds a very balanced tree
    void rebuild();
    // used after a rebuild to restore color bits
    void rb_recolor(bool);

    // red-black repair function used to restore rb-ness on insert
    void rb_repair();

    // an insert which will allow for the creation of a new bitdictionary.
    //  used by insert when new bitdictionaries are needed.
    void insert_rec(word_t p,bool b,word_t &s);

    // just something to check up on the tree structure
    void depth_maxmin(word_t &,word_t &) const;
  };

  // definition of dictionary, just an vnode which is
  //  responsible for its deletion
  template <class BitDict>
  class dictionary_T: public vnode_T<BitDict> {
  public:
    dictionary_T(): vnode_T<BitDict>() {
      vnode_T<BitDict>::_sum = 0;
      vnode_T<BitDict>::_size = 0;
      vnode_T<BitDict>::_is_leaf = 1;
      vnode_T<BitDict>::data.u_leaf = BitDict::new_bits();
    }
    // delete all allocations
    ~dictionary_T() { vnode_T<BitDict>::destroy(); }
    void clear() {
      vnode_T<BitDict>::destroy();
      vnode_T<BitDict>::_sum = 0;
      vnode_T<BitDict>::_size = 0;
      vnode_T<BitDict>::_is_leaf = 1;
      vnode_T<BitDict>::data.u_leaf = BitDict::new_bits();
    }
  }; // class dictionary
};

// node
template <class BitDict>
BW::vnode_T<BitDict>::node::node(typename BitDict::bits *l) {
  rb=0;
  L_is_leaf=1; R_is_leaf=1;
  L.u_leaf = l; R.u_leaf = BitDict::new_bits();

  word_t  smL;
  word_t szL;
  BitDict::divide(L.u_leaf,szL,R.u_leaf,smL);
  sumL  = smL;
  sizeL = szL;
}

// vnode methods
template <class BitDict>
void BW::vnode_T<BitDict>::destroy() {
  if (is_leaf()) {
    BitDict::delete_bits(data.u_leaf);
  }
  else {
    vnode_T<BitDict> l(*this),r(*this);
    l.goleft().destroy();
    r.goright().destroy();
    delete data.u_node;
  }
}

template <class BitDict>
word_t BW::vnode_T<BitDict>::find0(word_t r) const {// find nth 0
  vnode_T<BitDict> i(*this);
  word_t ret=0;
  while(!i.is_leaf()) {
    if (i.data.u_node->sizeL - i.data.u_node->sumL < r) {
      r -= i.data.u_node->sizeL - i.data.u_node->sumL;
      ret += i.data.u_node->sizeL;
      i.goright();
    }
    else {
      i.goleft();
    }
  }
  ret+=BitDict::find0(i.data.u_leaf,i.size(),r);
  return ret;
}

template <class BitDict>
word_t BW::vnode_T<BitDict>::find1(word_t r) const {// find nth 1
  vnode_T<BitDict> i(*this);
  word_t ret=0;
  while(!i.is_leaf()) {
    if (i.data.u_node->sumL < r) {
      r -= i.data.u_node->sumL;
      ret += i.data.u_node->sizeL;
      i.goright();
    }
    else {
      i.goleft();
    }
  }
  ret+=BitDict::find1(i.data.u_leaf,i.size(),r);
  return ret;
}

template <class BitDict>
void BW::vnode_T<BitDict>::lookup(word_t p,word_t &s,bool &b) const {
  vnode_T<BitDict> i(*this);
  s = 0;
  while(!i.is_leaf()) {
    if (p < i.data.u_node->sizeL) {
      i.goleft();
    }
    else {
      p -= i.data.u_node->sizeL;
      s += i.data.u_node->sumL;
      i.goright();
    }
  }
  BitDict::lookup(i.data.u_leaf,i.size(),p,s,b);
}

template <class BitDict>
void BW::vnode_T<BitDict>::lookup(word_t p,word_t &s) const {
  vnode_T<BitDict> i(*this);
  s = 0;
  while(!i.is_leaf()) {
    if (p < i.data.u_node->sizeL) {
      i.goleft();
    }
    else {
      p -= i.data.u_node->sizeL;
      s += i.data.u_node->sumL;
      i.goright();
    }
  }
  BitDict::lookup(i.data.u_leaf,i.size(),p,s);
}

template <class BitDict>
void BW::vnode_T<BitDict>::lookup(word_t p,bool &b) const {
  vnode_T<BitDict> i(*this);
  while(!i.is_leaf()) {
    if (p < i.data.u_node->sizeL) {
      i.goleft();
    }
    else {
      p -= i.data.u_node->sizeL;
      i.goright();
    }
  }
  BitDict::lookup(i.data.u_leaf,i.size(),p,b);
}

template <class BitDict>
void BW::vnode_T<BitDict>::insert(word_t p,bool b,word_t &s) {
  vnode_T<BitDict> i(*this);
  word_t savep=p;
  s = 0;

  // increment total counts in top level
  _size++; _sum += b;
  while(!i.is_leaf()) {
    if (p < i.data.u_node->sizeL) {
      i.data.u_node->sizeL++;
      i.data.u_node->sumL+=b;
      i.goleft();
    }
    else {
      p -= i.data.u_node->sizeL;
      s += i.data.u_node->sumL;
      i.goright();
    }
  }
  if (BitDict::insert(i.data.u_leaf,i.size(),p,b,s)) {
    return;
  }

  // repair total counts in top level
  _size--; _sum -= b;

  p = savep;
  s = 0;

  i = (*this);
  while(!i.is_leaf()) {
    if (p+1 < i.data.u_node->sizeL) {
      i.data.u_node->sizeL--;
      i.data.u_node->sumL-=b;
      i.goleft();
    }
    else if (p+1 > i.data.u_node->sizeL) {
      p -= i.data.u_node->sizeL;
      i.goright();
    }
    else {
      fprintf(stderr,"insert_repair:: corruption detected\n");
      exit(9);
    }
  }
  // now do the recursive insert
  insert_rec(savep,b,s);
}

template <class BitDict>
void BW::vnode_T<BitDict>::insert_rec(word_t p,bool b,word_t &s) { // sum and insert bit value
  if (is_leaf()) {
    if (!BitDict::insert(data.u_leaf,size(),p,b,s)) {
      //fprintf(stderr,"DIVIDE!!! %u %u\n",size(),sum());
      // oh no! we are full, time to repair
      node *neu = new node(data.u_leaf);

      //printf("new node %x inserted at %x\n",neu,data.u_leaf);
      // hookup to the parent
      if (parent) {
	if (parent->data.u_node->L.u_leaf == data.u_leaf) {
	  parent->data.u_node->L_is_leaf = 0;
	  parent->data.u_node->L.u_node = neu;
	}
	else {
	  parent->data.u_node->R_is_leaf = 0;
	  parent->data.u_node->R.u_node = neu;
	}
      }
      _is_leaf = 0;
      data.u_node = neu;

      // do the insert again on the two new leaves
      if (p < data.u_node->sizeL) {
	if (!BitDict::insert(data.u_node->L.u_leaf,data.u_node->sizeL,p,b,s))
	  {fprintf(stderr,"Weirdness 1\n"); exit(5);}
	data.u_node->sizeL++;
	data.u_node->sumL+=b;
      }
      else {
	p -= data.u_node->sizeL;
	s += data.u_node->sumL;
	if (!BitDict::insert(data.u_node->R.u_leaf,size()-data.u_node->sizeL,p,b,s))
	  {fprintf(stderr,"Weirdness 2\n"); exit(5);}
      }

      data.u_node->rb = 1;
      rb_repair();
    }
    else { // need an actual node
      // everything worked fine
    }
  }
  else {
    vnode_T<BitDict> i(*this);
    i.parent = this;
    if (p < data.u_node->sizeL) {
      i.goleft();
      data.u_node->sizeL++;
      data.u_node->sumL+=b;
      i.insert_rec(p,b,s);
    }
    else {
      i.goright();
      p -= data.u_node->sizeL;
      s += data.u_node->sumL;
      i.insert_rec(p,b,s);
    }
  }
  // increment total counts to be consistent
  _size++; _sum += b;
}

template <class BitDict>
void BW::vnode_T<BitDict>::rb_repair() {
  // rb=1 red, rb=0 black
  if (!parent) {
    // the root is always black
    data.u_node->rb = 0;
  }
  else if (data.u_node->rb && parent->data.u_node->rb) {
    vnode_T<BitDict> *grandparent = parent->parent;

    node *self= data.u_node;
    node *par = parent->data.u_node;
    node *gpar= grandparent->data.u_node;

    bool selfL = (self == par->L.u_node);
    bool parL  = (par  == gpar->L.u_node);

    // is my uncle red?
    if ((parL && !gpar->R_is_leaf && gpar->R.u_node->rb) ||
	(!parL && !gpar->L_is_leaf && gpar->L.u_node->rb) ) {
      //printf("doing case 1\n");
      gpar->rb = 1;
      gpar->R.u_node->rb = 0;      
      gpar->L.u_node->rb = 0;
      // the only place where we need a recursion
      grandparent->rb_repair();
    }
    else { // my uncle is black, we do rotation and balance
      if (parL) {
	//printf("case 2 rotation left\n");
	if (!selfL) {
	  par->R  = self->L;     par->R_is_leaf = self->L_is_leaf;
	  self->L.u_node = par;   self->L_is_leaf = 0;
	  gpar->L.u_node = self;
	  self->sizeL += par->sizeL;
	  self->sumL += par->sumL;
	  par = gpar->L.u_node;
	  self = par->L.u_node;
	}
	//printf("case 3 rotation left\n");
	par->L = par->R;  par->L_is_leaf = par->R_is_leaf;
	par->R = gpar->R; par->R_is_leaf = gpar->R_is_leaf;
	gpar->R.u_node = par;  gpar->R_is_leaf=0;
	gpar->L.u_node = self;

	word_t gpsize = gpar->sizeL;
	word_t  gpsum  = gpar->sumL;
	gpar->sizeL = par->sizeL;
	gpar->sumL  = par->sumL;
	par->sizeL = gpsize - par->sizeL;
	par->sumL  = gpsum - par->sumL;
      }
      else {
	//printf("case 3 rotation left\n");
	if (selfL) {
	  par->L  = self->R;   par->L_is_leaf = self->R_is_leaf;
	  self->R.u_node = par; self->R_is_leaf = 0;
	  gpar->R.u_node = self;
	  par->sizeL -= self->sizeL;
	  par->sumL -= self->sumL;
	  par = gpar->R.u_node;
	  self = par->R.u_node;
	}
	//printf("case 3 rotation right\n");
	par->R = par->L;  par->R_is_leaf = par->L_is_leaf;
	par->L = gpar->L; par->L_is_leaf = gpar->L_is_leaf;
	gpar->L.u_node = par;  gpar->L_is_leaf=0;
	gpar->R.u_node = self;

	word_t psize = par->sizeL;
	word_t  psum  = par->sumL;
	par->sizeL = gpar->sizeL;
	par->sumL  = gpar->sumL;
	gpar->sizeL += psize;
	gpar->sumL  += psum;
      }
    }
  }
}

template <class BitDict>
void BW::vnode_T<BitDict>::fdump(FILE *f) const {
  if (is_leaf()) {
    BitDict::fdump(data.u_leaf,f,size());
  }
  else {
    fprintf(f,"(%c" word_FMT ":" word_FMT ":",
	    char('q'+data.u_node->rb),data.u_node->sizeL,data.u_node->sumL);
    vnode_T<BitDict> l(*this),r(*this);
    fputc('L',f); l.goleft().fdump(f);
    fputc('R',f); r.goright().fdump(f);
    fprintf(f,")");
  }
};

template <class BitDict>
typename BW::vnode_T<BitDict>::node *BW::vnode_T<BitDict>::flatten(word_t &szR,
								   word_t &smR) {
  vnode_T<BitDict> le(*this),ri(*this);
  node  *endL,*endR;

  if (is_leaf()) {
    szR = size();
    smR = sum();

    return NULL;
    // flatten on a leaf does nothing, but gets the counts
  }
  
  // turn all flattened nodes red
  data.u_node->rb = 1;

  word_t szL;  word_t smL;
  endL = le.goleft().flatten(szL,smL);
  endR = ri.goright().flatten(szR,smR);

  if (endR) {
    data.u_node->R_is_leaf = 0;
    data.u_node->R         = ri.data;
  }
  else {
    endR = data.u_node;
  }
  if (endL) {
    data.u_node->sizeL = szL;
    data.u_node->sumL  = smL;

    data.u_node->L_is_leaf = endL->R_is_leaf;
    data.u_node->L         = endL->R;

    endL->R_is_leaf = 0;
    endL->R         = data;
  
    data.u_node = le.data.u_node;
  }
  return endR;
}

template <class BitDict>
void BW::vnode_T<BitDict>::rebuild() {
  if (is_leaf()) return;

  word_t nnodes=0;
  for(vnode_T<BitDict> i=*this; !i.is_leaf(); i.goright()) ++nnodes;

  word_t pow2,lg;
  for(pow2=1,lg=1; 2*pow2 <= nnodes; pow2 *= 2,++lg);

  for(word_t lf=nnodes - pow2; pow2>1 ; pow2/=2, lf = pow2-1) {
    vnode_T<BitDict> i = (*this);
    for(word_t t=0; t < lf; ++t) {
      // leapfrog move
      // i.data becomes new "root" and i+1.data its child
      node *id = i.data.u_node;
      node *i1 = id->R.u_node;
      word_t  sumi1 = i1->sumL;
      word_t sizei1= i1->sizeL;
      
      // i leapfrogs i1
      id->R.u_node = i1->R.u_node;
      // i1 passes L to R and takes i's L
      i1->R_is_leaf = i1->L_is_leaf;
      i1->R         = i1->L;
      i1->L_is_leaf = id->L_is_leaf;
      i1->L         = id->L;
      
      // i's L now points to i1
      id->L_is_leaf = 0;
      id->L.u_node  = i1;
      
      // adjust counts
      i1->sizeL   = id->sizeL;
      i1->sumL    = id->sumL;
      id->sizeL   += sizei1;
      id->sumL    += sumi1;
      
      // now just go right...
      i.goright();
    }
  }

  //word_t mx=0,mn=100;
  //depth_maxmin(mx,mn);
  //fprintf(stderr,"depths = %lu %lu\n",mn,mx);
  rb_recolor( (lg % 2) );
  this->data.u_node->rb=0; // root is always black
}

template <class BitDict>
void BW::vnode_T<BitDict>::depth_maxmin(word_t &mx,word_t &mn) const {
  if (is_leaf()) {
    mx = 1; mn = 1;
  }
  else {
    word_t m1,n1,m2,n2;
    vnode_T<BitDict> L(*this),R(*this);
    L.goleft().depth_maxmin(m1,n1);
    R.goright().depth_maxmin(m2,n2);
    mx = ((m1>m2)?m1:m2)+1;
    mn = ((n1<n2)?n1:n2)+1;
  }
};

template <class BitDict>
void BW::vnode_T<BitDict>::rb_recolor(bool col) {
  if (is_leaf()) {
    ;
  }
  else {
    data.u_node->rb = col;
    vnode_T<BitDict> L(*this),R(*this);
    L.goleft().rb_recolor(!col);
    R.goright().rb_recolor(!col);
  }
}

template <class BitDict>
bool BW::vnode_T<BitDict>::regrow(bool lastrb) {
  if (is_leaf()) {
    return true;
  }
  else {
    if (lastrb && data.u_node->rb) {
      rb_repair();
      return false;
    }
    else {
      vnode_T<BitDict> ri(*this);
      ri.parent = this;
      return ri.goright().regrow(data.u_node->rb);
    }
  }
}

template <class BitDict>
void BW::vnode_T<BitDict>::merge() {
  // looks icky, but we're trying to avoid stack growth here
  if (!is_leaf()) {
    node *oldn=NULL,*n=data.u_node;
    word_t rsize=size();
    word_t  rsum =sum();

    while(!n->R_is_leaf) {// eat from L children of R
      node *m = n->R.u_node;
      if (!m->sizeL) {
	// leapfrog over m
	n->R_is_leaf = m->R_is_leaf; n->R = m->R;
	// destroy m and its L child
	BitDict::delete_bits(m->L.u_leaf);
	delete m;
      }
      else {
	word_t  nszL=n->sizeL,mszL=m->sizeL;
	word_t   nsmL=n->sumL,msmL=m->sumL;
	if (BitDict::merge(n->L.u_leaf, nszL, nsmL,
				 m->L.u_leaf, mszL, msmL  )) {
	  n->sizeL=nszL; n->sumL =nsmL;
	  m->sizeL=mszL; m->sumL =msmL;
	}
	else {
	  // can't eat anymore, so move right
	  rsize-= nszL; rsum-= nsmL; 
	  oldn=n;  n=m;
	}
      }
    }

    // finally eat from the very last leaf
    // compute the leaf's  size and sum
    word_t  nszL=n->sizeL; word_t   nsmL=n->sumL;
    rsize -= nszL; rsum  -= nsmL;
  
    BitDict::merge(n->L.u_leaf, nszL, nsmL,
			 n->R.u_leaf, rsize, rsum );

    n->sizeL=nszL; n->sumL =nsmL;
    if (!rsize) {
      BitDict::delete_bits(n->R.u_leaf);
      if (oldn) {
	oldn->R_is_leaf = 1;
	oldn->R.u_leaf = n->L.u_leaf;
      }
      else { // I am in the root, so I should patch that up
	_is_leaf = 1;
	data.u_leaf = n->L.u_leaf;
      }
      delete n;
    }
  }
}

template <class BitDict>
void BW::vnode_T<BitDict>::freeze() {
  word_t s;  word_t sm;
  flatten(s,sm);
  merge();
  rebuild();
  //while(!regrow(1)) ;
}

template <class BitDict>
void BW::vnode_T<BitDict>::thaw() {
  word_t bh=0;
  if (!rb_check(bh)) {
    fprintf(stderr,"WARNING: red-black property being restored\n");
    word_t s=0;  word_t sm=0;
    flatten(s,sm);
    while(!regrow(1)) ;
  }
}

// useful for testing rb tree stuff
template <class BitDict>
word_t BW::vnode_T<BitDict>::depth() const {
  if (is_leaf()) {
    return 0;
  }
  else {
    vnode_T<BitDict> l(*this),r(*this);
    word_t dl=1+l.goleft().depth();
    word_t dr=1+r.goright().depth();
    if (dr > dl) return dr;
    else return dl;
  }
}
template <class BitDict>
word_t BW::vnode_T<BitDict>::num_nodes() const {
  if (is_leaf()) {
    return 0;
  }
  else {
    vnode_T<BitDict> l(*this),r(*this);
    return 1+l.goleft().num_nodes()+r.goright().num_nodes();
  }
}

template <class BitDict>
bool BW::vnode_T<BitDict>::rb_check(word_t &bh) const {
  if (is_leaf()) {
    bh = 0;
    return true;
  }
  else {
    vnode_T<BitDict> l(*this),r(*this);
    word_t dl=0,dr=0;
    if (l.goleft().rb_check(dl) && r.goright().rb_check(dr)) {
      if (data.u_node->rb &&
	  (!l.is_leaf() && l.data.u_node->rb)) 
	{fprintf(stderr,"red node has red left child\n"); return false;}
      if (data.u_node->rb &&
	  (!r.is_leaf() && r.data.u_node->rb))
	{fprintf(stderr,"red node has red right child\n"); return false;}
      if (dl == dr) {
	bh = dl;
	bh += (1-data.u_node->rb);
	return true;
      }
      else {
	{fprintf(stderr,"inconsistent black depth\n"); return false;}
      }
    }
    else {
      return false;
    }
  }
}

template <class BitDict>
bool BW::vnode_T<BitDict>::fprint(bitio &bf) const {
  if (is_leaf()) {
    return BitDict::fprint(data.u_leaf,bf,size());
  }
  else {
    vnode_T<BitDict> l(*this),r(*this);
    if (!l.goleft().fprint(bf)) return false;
    if (!r.goright().fprint(bf)) return false;
  }
  return true;
}

template <class BitDict>
bool BW::vnode_T<BitDict>::fprint(FILE *f) const {
  word_t s = size();
  word_t  sm= sum();

  if (!fwrite(&s,   sizeof(word_t),    1,f)) return false;
  if (!fwrite(&sm,   sizeof(word_t),    1,f)) return false;

  bitio bf(f,"w");

  if (!fprint(bf)) return false;

  return bf.flush();
}

template <class BitDict>
bool BW::vnode_T<BitDict>::fscan(bitio &bf) {
  node *par = NULL;
  vnode_T<BitDict> ptr(*this);

  if (!is_leaf()) {
    fprintf(stderr,"BWdictionary: unable to fscan if dictionary not clear()-ed\n");
    exit(9);
  }

  for(;;) {
    word_t sz = ptr.size();
    word_t  sm = ptr.sum();
    if (!BitDict::fscan(ptr.data.u_leaf,bf,sz,sm)) return false;

    if (sz != ptr.size()) { // crap! I need to be a node and pass along
      node *neu = new node();

      neu->rb = 1;
      neu->sizeL = sz;
      neu->sumL = sm;
      neu->L_is_leaf = 1;
      neu->L.u_leaf = ptr.data.u_leaf;
      neu->R_is_leaf = 1;
      neu->R.u_leaf = BitDict::new_bits();

      if (par) {
	par->R_is_leaf = 0;
	par->R.u_node =  neu;
      }
      else {
	this->_is_leaf = 0;
	this->data.u_node = neu;
      }

      ptr._is_leaf = 0;
      ptr.data.u_node = neu;
      ptr.goright();
      
      par = neu;
    }
    else { // I got just what I needed, whew!
      if (sm != ptr.sum()) {fprintf(stderr,"Weirdness 7: checksum failed\n"); exit(7);}
      break;
    }
  } // for

  rebuild(); // produces a better depth tree than regrow
  //while(!regrow(1)) ;
  return true;
}

template <class BitDict>
bool BW::vnode_T<BitDict>::fscan(FILE *f) {
  word_t s ;
  word_t  sm;
  if (!fread(&s,   sizeof(word_t),    1,f)) return false;
  if (!fread(&sm,  sizeof(word_t),     1,f)) return false;

  _size = s;
  _sum  = sm;
  bitio bf(f,"r");

  return fscan(bf);
}


#endif

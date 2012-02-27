/***************************************************************************
 * Title:          PagedList.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PAGED_LIST_H_
#define PAGED_LIST_H_

#include <iostream>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include "utils.h"
using namespace std;
// An implementation of a paged linked list. 
// Each page stores 'SIZE' elemnets, and pages are linked in a list.
// This allows fast append operations, and doesn't have to resize memory.


template <typename T, ssize_t SIZE> 
class Page {
public:
	ssize_t  size;
  T page[SIZE];
  Page *next;
	//  Page *head;
	Page() {
		size = 0;
		next = NULL;
	}
};

template <typename T, ssize_t SIZE>
class myterator {
public:
  Page<T,SIZE> *list;
  ssize_t pos;
  myterator() {
    list = NULL;
    pos = 0;
  }
  myterator(Page<T,SIZE> *l, ssize_t p) {
    list = l; pos = p;
  }

  bool operator==(const myterator &it) const {
    return (it.list == list and pos == it.pos);
  }
	
	bool operator!=(const myterator &it) const {
		return it.list != list or pos != it.pos;
	}

  myterator & operator=(const myterator &it) {
		if (this != &it) {
			list = it.list;
			pos  = it.pos;
		}
    return *this;
  }
  T& operator*() {
    assert(list != NULL);
    assert(pos < SIZE);
    return (*list).page[pos];
  }
  myterator & operator++() {
    pos++;
    if (pos == SIZE) {
      list = list->next;
      pos  = 0;
    }
    return *this;
  }
  ssize_t AtEndOfPage() {
    return (list==NULL or pos == SIZE-1);
  }
  bool operator<(const myterator &it) const {
    if (it.list == NULL) 
      return 0;
    else if (it.list == list) {
      return pos < it.pos;
    }
    else {
      // it is on a different page than this one if this one is fist from 'head'
      // return 1
      Page<T,SIZE> *page;
      page = list->head;
      while (page != NULL and page != list and page != it.list) {
				page = page->next;
      }
      if (page == list) 
				return 1;
      else if (page == it.list)
				return 0;
      else {
				std::cout << "INTERNAL ERROR, Should have found a page" << std::endl;
				exit(1);
      }
    }
  }
};


template <typename T, ssize_t SIZE>
class PagedList {
public:
  Page<T,SIZE> *head, *next;
  Page<T,SIZE> *cur;
	T* curPtr;
  ssize_t lastPos;
  // Define an myterator class for packaged traversal of this 
  // data structure.  This is called an myterator, and looks like the stdlib c++ 
  // myterators, but it can't be used as one.  (maybe I should call it something 
  // different... myterator I think.
  typedef myterator<T, SIZE> iterator;
  T* End() {
    //iterator it(cur, last.pos);
		return &cur.page[lastPos];
		//    return it;
  }
  T* Begin() {
		//    iterator it(head, 0);
		if (head == NULL)
			return NULL;
		else
			return &(head->page)[0];
  }

private:
  // We don't want anybody to store anywhere
  iterator last;
  T* &Store(T &value) {
		*curPtr = value;

		// Update the n elements in the list,
		// and advance place of where to store the next element.
		++curPtr;
		++cur->size;
		//		cout << "sorting..." << endl;
		std::sort((*cur).page, &((*cur).page[cur->size]));		
		return curPtr;
  }
public:
  ssize_t numAdded;
  ssize_t numPages;
  PagedList() {
    head = cur = NULL;
    last.list = NULL; last.pos = 0;
		//UNUSED// T* curPtr = NULL;
    lastPos = 0;
    numAdded = 0;
    numPages = 0;
  }
  T* Find(T &value) {
		/*
			iterator it;
    it.list = head;
		it.pos  = 0;
    iterator end;
		end.list = cur;
		end.pos  = last.pos;
		*/
		Page<T,SIZE> *page = head;
		//UNUSED// ssize_t found = 0;
		T* location = NULL;
		T* lb;
		while(page != NULL) {
			//			cout << "searching: " << value << endl;
			//lb = std::lower_bound( &(*page).page[0], &(*page).page[page->size], value);
			lb = std::find( &(*page).page[0], &(*page).page[page->size], value);

			/*

			ssize_t i, ps  = page->size;
			T val;
			T* valPtr = &page->page[0], *pageEnd = &(*page).page[page->size];

				for (;valPtr != pageEnd and value < *valPtr; valPtr++);
							
				if (valPtr != pageEnd and value == *valPtr) {
				//				cout << "found: " << value << endl;
				location = valPtr;
				break;
			}
			*/
			if (lb != &(page->page[page->size]) and
					*lb == value) {
				location = lb;
				break;
			}
			else {
				page = page->next;
			}
		}
		return location;
  }

  ssize_t size() const {
    return numAdded;
  }
  double load() const {
    return double(numAdded) / numPages;
  }

	void Free() {
		Page<T,SIZE> *fcur, *fnext;
		fcur = head;
		if (fcur != NULL)
			next = fcur->next;
		while (fcur != NULL) {
			delete fcur;
			fcur  = fnext;
			fnext = fcur->next;
		}
		head = next = cur = NULL;
	}

  T* Append(T &value) {
    // Expand the list if necessary
		if (head == NULL) {
      Page<T,SIZE> *newPage = new Page<T,SIZE>;
			cur = head = newPage;
			curPtr = &(head->page[0]);
		}
    else if (curPtr ==  &((*cur).page[SIZE]) ) {
			// Stored the full list.  Need to process
			// this list, then create a new list to store into.
			/*			if (cur != NULL) {
				std::sort((*cur).page, &((*cur).page[SIZE]));
			}
			*/
      Page<T,SIZE> *newPage = new Page<T,SIZE>;
      ++numPages;
			cur->next = newPage;

			// Advance ptr to new list.
      cur = newPage;
			curPtr = &(*cur).page[0];

    }

    assert(head != NULL);
    // Store at the current last position
		//   assert(last.list != NULL);
    // Advance the last position
		//		cout << "about to store " << value << endl;
		Store(value);
    return curPtr;
  }
	
	T* Advance() {
		// set the pointer to the next position where it is ok.
		++curPtr;
	}
  friend std::ostream &operator<<(std::ostream &out, const PagedList &p) {
		//    iterator it, end;
		//    end = p.End();
    out << p.size() << std::endl;
		Page<T,SIZE> *curPage = p.head;
		while (curPage != NULL) {
			T* pageEnd = &(*curPage).page[curPage->size];
			T* it;
			for (it = &(*curPage).page[0]; it != pageEnd; ++it) {
				out << *it << std::endl;
			}
			curPage = curPage->next;
    }
    return out;
  }
};

    
#endif

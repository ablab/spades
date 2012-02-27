/***************************************************************************
 * Title:          PagedHashTable.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PAGED_HASH_TABLE
#define PAGED_HASH_TABLE

#include <vector>
#include "PagedList.h"


template<typename T, ssize_t PAGE_SIZE, typename HashFunct, typename UpdateFunct>
class PagedHashTable {
public:
  ssize_t size;
	ssize_t count;
  PagedList<T,PAGE_SIZE> *table;
  typedef typename PagedList<T,PAGE_SIZE>::iterator iterator;
	typedef Page<T,PAGE_SIZE> HashPage;
	typedef PagedList<T,PAGE_SIZE> HashPagedList;
	typedef T Data;
	PagedHashTable() { count = 0;};
  PagedHashTable(HashFunct &hashFunct) {
		count = 0;
    size = hashFunct.MaxHashValue();
    table = new HashPagedList[size];
  }
	void Init(HashFunct &hashFunct) {
		size = hashFunct.MaxHashValue();
		table = new HashPagedList[size];
		//UNUSED// ssize_t i;
	}

  ssize_t CountSize() {
    ssize_t i;
    ssize_t total = 0;
    for (i = 0; i < size; i++ ) 
      total += table[i].numAdded;

    return total;
  }
  ssize_t CountPages() {
    ssize_t i;
    ssize_t numPages = 0;
    for (i = 0; i < size; i++ ) {
      numPages += table[i].numPages;
    }
    return numPages;
  }

  ssize_t Store(T value, HashFunct hashFunct = HashFunct(), UpdateFunct update = UpdateFunct()) {
    ssize_t index = hashFunct(value);
    if (index < 0) 
     return 0;
		//    typename PagedList<T,PAGE_SIZE>::iterator it;
		T* it;
    it = table[index].Find(value);
		//    if (it == table[index].End()) {
		if (it == NULL) {
			//			std::cout << "appending value: " << value << std::endl;
      table[index].Append(value);
			count++;
      return 1;
    }
    else {
      update(*it);
      return 0;
    }
  }
	
	ssize_t Find(T &value, HashFunct hashFunct=  HashFunct()) {
		ssize_t index = hashFunct(value);
		//    typename PagedList::iterator it;
		T* it;
		it = table[index].Find(value);
		//		if ( it == table[index].End()) {
		if (it == NULL) {
			// Store auxilliary data into value
			//			std::cout << "did not find " << value << std::endl;
			return 0;
		}
		else {
			//			std::cout << "found value: " << std::endl;
			value = *it;
			return 1;
		}
	}
		
	void Summarize() {
		ssize_t nChains = 0;
		ssize_t chainLength, totalChainLength = 0;
		ssize_t i;
		ssize_t maxLen = 50;
		ssize_t bins[maxLen];
		ssize_t b;
		for (b = 0; b < maxLen; b++ ) bins[b] = 0;
		HashPage *page;
		for (i =0 ; i < size; i++) {
			if (table[i].head != NULL) {
				nChains++;
				chainLength = 1;
				page = table[i].head;
				while (page != NULL) {
					page = page->next;
					chainLength++;
					if (chainLength  > maxLen-1) 
						bins[maxLen-1]++;
					else 
						bins[chainLength]++;
				}
				totalChainLength += chainLength;
			}
		}
		for (b = 0; b < maxLen; b++ )
			std::cout << bins[b] << " ";
		std::cout << nChains << " chains." << std::endl;
		std::cout << totalChainLength / nChains << " average chain length. " << std::endl;
	}

	void Free() {
		ssize_t i;
		for (i = 0; i < size; i++) {
			table[i].Free();
		}
	}
};

/*
  To finish later if I ever need it.

template<typename T, ssize_t PAGE_SIZE>
class PagedHashTableIterator {
  ssize_t index;
  typename PagedList<T,PAGE_SIZE>::iterator it;
public:
  T & operator*() {
  }
  bool operator==(const PagedHashTableItertor<T, PAGE_SIZE> &comp) const {
    return (it.index == comp.index and
	    it.list  == comp.list and
	    it.pos   == comp.pos);
  }
  typename PagedHashTableIterator & operator++() {
    if (it == End())
      return it;
  }
  
};
*/

    
    

#endif

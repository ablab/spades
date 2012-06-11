/*
** TrieFactory.h
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:02:33 2006 Julien Lemoine
** $Id$
** 
** Copyright (C) 2006 Julien Lemoine
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
** 
** You should have received a copy of the GNU Lesser General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef   	TRIEFACTORY_H_
# define   	TRIEFACTORY_H_

#include <list>

namespace ToolBox
{
  // fwd declaration
  template <typename T>
  class TrieNode;

  /**
   * The goal of this class is to allocate Trie node by paquet of X
   * element in order to reduce heap-admin size
   */
  template <typename T>
  class TrieFactory
    {
    public:
      TrieFactory(unsigned paquetSize);
      ~TrieFactory();

    private:
      /// avoid default constructor
      TrieFactory();
      /// avoid copy constructor
      TrieFactory(const TrieFactory &e);
      /// avoid affectation operator
      TrieFactory& operator=(const TrieFactory &e);

    public:
      TrieNode<T>* getNewNode(const T &value);
      void clear();

    private:
      unsigned			_paquetSize;
      std::list<TrieNode<T>*>	_allocatedNodes;
      TrieNode<T>		*_lastNodes;
      unsigned			_nbUsedInLastNodes;
    };
}

#endif 	    /* !TRIEFACTORY_H_ */

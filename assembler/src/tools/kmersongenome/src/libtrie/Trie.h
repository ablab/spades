/*
** Trie.h
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:42:53 2006 Julien Lemoine
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

#ifndef   	TRIE_H_
# define   	TRIE_H_

#include <iostream>

namespace ToolBox
{
  //fwd declaration
  template <typename T>
    class TrieFactory;
  template <typename T>
    class TrieNode;

  /**
   * Implement a trie in memory with the smallest structure as possible
   * (use few RAM as possible)
   */
  template <typename T>
    class Trie
    {
    public:
      /// constuctor, empty is the value returned when no match in found
      /// in trie
      Trie(const T &empty);
      ~Trie();

    private:
      /// avoid default constructor
      Trie();
      /// avoid copy constructor
      Trie(const Trie &e);
      /// avoid affectation operator
      Trie& operator=(const Trie &e);
    
    public:
      /// add an entry in the Trie, if entry already exist an exception
      /// is throw
      void addEntry(const char *str, unsigned strLen, const T &value);
      /// associates a value to a string, if string is already in Trie,
      /// value is overwriten
      void setEntry(const char *str, unsigned strLen, const T &value);
      /// get an entry in the Trie
      const T& getEntry(const char *str, unsigned strLen) const;
      /// get initial TrieNode
      const TrieNode<T>* getInitialNode() const;
      /// display content of trie in output stream
      void display(std::ostream &os);
      /// clear the content of trie
      void clear();

    protected:
      TrieNode<T>* _addEntry(const char *str, unsigned strLen);

    private:
      /// value returned when no match is found in trie
      T			_empty;
      /// factory
      TrieFactory<T>	*_factory;
      /// first node of trie
      TrieNode<T>		*_initialNode;
    };
}

#endif 	    /* !TRIE_H_ */

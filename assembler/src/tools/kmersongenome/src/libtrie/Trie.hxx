/*
** Trie.hxx
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:47:20 2006 Julien Lemoine
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

#ifndef   	TRIE_HXX_
# define   	TRIE_HXX_

#include "Trie.h"
#include <vector>
#include <algorithm>
#include <assert.h>
#include "TrieFactory.hxx"
#include "TrieNode.hxx"
#include "Exception.h"

template <typename T>
ToolBox::Trie<T>::Trie(const T &empty) :
  _empty(empty), _factory(0x0), _initialNode(0x0)
{
  // initialize nodes by paquets of 10000
  _factory = new TrieFactory<T>(10000);
  _initialNode = _factory->getNewNode(_empty);
}

template <typename T>
ToolBox::Trie<T>::~Trie()
{
  if (_factory)
    delete _factory;
}

template <typename T>
void ToolBox::Trie<T>::setEntry(const char *str, unsigned strLen, const T &value)
{
  TrieNode<T>	*node = _addEntry(str, strLen);
  node->setValue(value);
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::Trie<T>::_addEntry(const char *str, unsigned strLen)
{
  unsigned	pos = 0;
  bool		found = true;
  TrieNode<T>	*node = _initialNode, *previous = 0x0;

  // Look for the part of the word which is in Trie
  while (found && pos < strLen)
    {
      found = false;
      previous = node;
      node = node->getSubNodeByLabel(str[pos]);
      if (node)
	{
	  found = true;
	  ++pos;
	}
    }

  // Add part of the word which is not in Trie
  if (!node || pos != strLen)
    {
      node = previous;
      for (unsigned i = pos; i < strLen; ++i)
	{
	  TrieNode<T> *newNode = _factory->getNewNode(_empty);
	  node->addSubNode(str[pos], newNode);
	  node = newNode;
	  ++pos;
	}
    }
  assert(node != 0x0);
  return node;
}

template <typename T>
void ToolBox::Trie<T>::addEntry(const char *str, unsigned strLen, const T &value)
{
  TrieNode<T>	*node = _addEntry(str, strLen);

  // Set the value on the last node
  if (node && node->getValue() != _empty)
    throw ToolBox::Exception("The word is already in automaton", HERE);
  node->setValue(value);
}

template <typename T>
const T& ToolBox::Trie<T>::getEntry(const char *str, unsigned strLen) const
{
  unsigned		pos = 0;
  bool			found = true;
  const TrieNode<T>	*node = _initialNode;
	
  while (found && pos < strLen)
    {
      found = false;
      node = node->getSubNodeByLabel(str[pos]);
      if (node)
	{
	  found = true;
	  ++pos;
	}
    }
  if (node && pos == strLen) // The word is complet in the automaton
    return node->getValue();
  return _empty;
}

template <typename T>
const ToolBox::TrieNode<T>* ToolBox::Trie<T>::getInitialNode() const
{
  return _initialNode;
}

template <typename T>
void ToolBox::Trie<T>::clear()
{
  _factory->clear();
  _initialNode = _factory->getNewNode(_empty);
}

template <typename T>
void ToolBox::Trie<T>::display(std::ostream &os)
{
  if (_initialNode)
    _initialNode->display(os, 0, 0);
}

#endif	    /* !TRIE_HXX_ */

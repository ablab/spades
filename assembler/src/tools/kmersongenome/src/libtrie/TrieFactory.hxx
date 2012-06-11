/*
** TrieFactory.hxx
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:07:49 2006 Julien Lemoine
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

#ifndef   	TRIEFACTORY_HXX_
# define   	TRIEFACTORY_HXX_

#include "TrieFactory.h"
#include "TrieNode.hxx"

template <typename T>
ToolBox::TrieFactory<T>::TrieFactory(unsigned paquetSize) :
  _paquetSize(paquetSize), _lastNodes(0x0), _nbUsedInLastNodes(0)
{
  _lastNodes = new TrieNode<T>[paquetSize];
}

template <typename T>
ToolBox::TrieFactory<T>::~TrieFactory()
{
  typename std::list<TrieNode<T>*>::const_iterator it;

  for (it = _allocatedNodes.begin(); it != _allocatedNodes.end(); ++it)
    delete[] *it;
  if (_lastNodes)
    delete[] _lastNodes;
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::TrieFactory<T>::getNewNode(const T &value)
{
  if (_nbUsedInLastNodes == _paquetSize)
    {
      _allocatedNodes.push_back(_lastNodes);
      _nbUsedInLastNodes = 0;
      _lastNodes = new TrieNode<T>[_paquetSize];
    }
  TrieNode<T> *res = &_lastNodes[_nbUsedInLastNodes];
  ++_nbUsedInLastNodes;
  res->setValue(value);
  res->clear();
  return res;
}

template <typename T>
void ToolBox::TrieFactory<T>::clear()
{
  typename std::list<TrieNode<T>*>::const_iterator it;
  for (it = _allocatedNodes.begin(); it != _allocatedNodes.end(); ++it)
    delete[] *it;
  _allocatedNodes.clear();
  _nbUsedInLastNodes = 0;
}

#endif	    /* !TRIEFACTORY_HXX_ */

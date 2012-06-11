/*
** TrieNode.hxx
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:23:30 2006 Julien Lemoine
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

#ifndef   	TRIENODE_HXX_
# define   	TRIENODE_HXX_

#include "TrieNode.h"

template <typename T>
ToolBox::TrieNode<T>::TrieNode() :
  _brother(0), _brotherLabel(0), _firstSubNode(0), _firstSubNodeLabel(0)
  /// we can not set _value here because type is unknown. assert that
  /// the value is set later with setValue()
{
}

template <typename T>
ToolBox::TrieNode<T>::~TrieNode()
{
  // do not delete _brother and _firstSubNode because they are
  // allocated by factory (TrieFactory) and factory will delete them
}

template <typename T>
void ToolBox::TrieNode<T>::setValue(const T &val)
{
  _value = val;
}

template <typename T>
const T& ToolBox::TrieNode<T>::getValue() const
{
  return _value;
}

template <typename T>
const ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getBrother() const
{
  return _brother;
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getBrother()
{
  return _brother;
}

template <typename T>
const ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::_getBrother(unsigned char chr) const
{
  const TrieNode<T> *brother = _brother;
  return _sequentialSearch(brother, _brotherLabel, chr);
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::_getBrother(unsigned char chr)
{
  return _sequentialSearch(_brother, _brotherLabel, chr);
}

template <typename T>
void ToolBox::TrieNode<T>::_addBrother(unsigned char chr, TrieNode<T> *brother)
{
  if (!_brother || _brotherLabel > chr)
    {
      brother->_setBrother(_brother, _brotherLabel);
      _brother = brother;
      _brotherLabel = chr;
    }
  else
    _brother->_addBrother(chr, brother);
}

template <typename T>
unsigned char ToolBox::TrieNode<T>::getBrotherLabel() const
{
  return _brotherLabel;
}

template <typename T>
const ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getSubNode() const
{
  return _firstSubNode;
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getSubNode()
{
  return _firstSubNode;
}

template <typename T>
unsigned char ToolBox::TrieNode<T>::getSubNodeLabel() const
{
  return _firstSubNodeLabel;
}

template <typename T>
const ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getSubNodeByLabel(unsigned char chr) const
{
  const TrieNode<T> *first = _firstSubNode;
  return _sequentialSearch(first, _firstSubNodeLabel, chr);
}

template <typename T>
ToolBox::TrieNode<T>* ToolBox::TrieNode<T>::getSubNodeByLabel(unsigned char chr)
{
  return _sequentialSearch(_firstSubNode, _firstSubNodeLabel, chr);
}

template <typename T>
void ToolBox::TrieNode<T>::addSubNode(unsigned char chr, TrieNode<T> *node)
{
  if (!_firstSubNode || _firstSubNodeLabel > chr)
    {
      node->_setBrother(_firstSubNode, _firstSubNodeLabel);
      _firstSubNode = node;
      _firstSubNodeLabel = chr;
    }
  else
    _firstSubNode->_addBrother(chr, node);
}

template <typename T>
template <typename Node>
inline Node ToolBox::TrieNode<T>::_sequentialSearch(Node first, unsigned char label, unsigned char val) const
{
  if (first && label <= val)
    {
      if (label == val)
	return first;
      return first->_getBrother(val);
    }
  return 0x0;
}

template <typename T>
void ToolBox::TrieNode<T>::_setBrother(TrieNode<T> *brother, unsigned char brotherLabel)
{
  _brother = brother;
  _brotherLabel = brotherLabel;
}

template <typename T>
void ToolBox::TrieNode<T>::display(std::ostream &os, unsigned offset, unsigned char label) const
{
  unsigned int i;
  for (i = 0; i < offset; ++i)
    os << " ";
  if (label)
    os << "label[" << label << "] ";
  os << "value[" << _value << "]" << std::endl;
  if (_firstSubNode)
    _firstSubNode->display(os, offset + 2, _firstSubNodeLabel);
  if (_brother)
    _brother->display(os, offset, _brotherLabel);
}

template <typename T>
void ToolBox::TrieNode<T>::clear()
{
  _brother = 0x0;
  _brotherLabel = 0;
  _firstSubNode = 0x0;
  _firstSubNodeLabel = 0;
}

#endif	    /* !TRIENODE_HXX_ */

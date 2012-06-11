/*
** TrieNode.h
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Sun Oct 29 16:15:57 2006 Julien Lemoine
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

#ifndef   	TRIENODE_H_
# define   	TRIENODE_H_

namespace ToolBox
{
  /**
   * @brief this class represent the node of a trie, it contains a
   * link to a sub node and a link to a brother (node which have the
   * same father)
   */
  template <typename T>
  class TrieNode
  {
  public:
    TrieNode();
    ~TrieNode();

  private:
    /// avoid copy constructor
    TrieNode(const TrieNode &e);
    /// avoid affectation operator
    TrieNode& operator=(const TrieNode &e);

  public:
    /// set value associed to node
    void setValue(const T &val);
    /// get value associed to node
    const T& getValue() const;

    /// get brother (return 0x0 this node has no brother)
    const TrieNode<T>* getBrother() const;
    TrieNode<T>* getBrother();
    /// get brother label
    unsigned char getBrotherLabel() const;

    // get first sub Node
    const TrieNode<T>* getSubNode() const;
    TrieNode<T>* getSubNode();
    unsigned char getSubNodeLabel() const;

    // Looking for a sub node
    const TrieNode<T>* getSubNodeByLabel(unsigned char chr) const;
    TrieNode<T>* getSubNodeByLabel(unsigned char chr);

    // add an edge
    void addSubNode(unsigned char chr, TrieNode<T> *node);

    /// display content of node in output stream
    void display(std::ostream &os, unsigned offset, unsigned char label) const;

    /// clear content of TrieNode
    void clear();

  protected:
    template <typename Node>
    Node _sequentialSearch(Node first, unsigned char label,
			   unsigned char val) const;
    /// set brother (used by sort)
    void _setBrother(TrieNode<T> *brother, unsigned char brotherLabel);
    /// add a new brother
    void _addBrother(unsigned char chr, TrieNode<T> *brother);
    /**
     * @ brief get brother that has the label chr (return 0x0 if brother is
     * not found)
     */
    const TrieNode<T>* _getBrother(unsigned char chr) const;
    /**
     * @ brief get brother that has the label chr (return 0x0 if brother is
     * not found)
     */
    TrieNode<T>* _getBrother(unsigned char chr);

  private:
    /// pointer to brother (node with same father as this one)
    TrieNode<T>		*_brother;
    /// character to go to brother (node with same father as this one)
    unsigned char	_brotherLabel;
    /// pointer to first sub node
    TrieNode<T>		*_firstSubNode;
    /// character to go to first subnode
    unsigned char	_firstSubNodeLabel;
    /// value associed to this node
    T			_value;
  };
}

#endif 	    /* !TRIENODE_H_ */

/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        tree.h                                                 *
* Source:      GNU C++                                                *
* Purpose:     tree headers                                           *
* History:     080928 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3TREEH__
#define __EV3TREEH__

#define RCS7 "$Id: tree.cxx,v 1.13 2005/09/08 23:37:23 liberti Exp liberti $"

#include<iostream>
#include<vector>
#include<string>
#include "common.h"
#include "exceptions.h"

template<class NodeType> class Pointer {

private:

  // pointer to node data
  NodeType* node;

  // how many Pointers point to the NodeType pointed to by *node
  Int* ncount;

public:

  // constructors  
  Pointer();
  Pointer(NodeType& n);
  Pointer(double v);
  Pointer(double c, int vi, std::string vn);
  Pointer(bool notinitialized);
  // copy constructor
  void SetTo(const Pointer<NodeType>& t);
  Pointer(const Pointer<NodeType>& t);
  // copy assignment: does NOT copy, it just references.
  Pointer<NodeType>& operator = (const Pointer<NodeType>& t);
  // copy factory
  // 1. this is a copy of pointer
  void SetToCopyOf(const Pointer<NodeType>& t);
  // 2. returns a copy of this
  Pointer<NodeType> Copy(void) const;
  // destructor
  void Destroy(void);
  ~Pointer();
  NodeType* operator->() const;
  NodeType GetPointee(void) const;
  bool operator == (const Pointer<NodeType>& t);
#ifdef MEMDEBUG
  NodeType* MemDebugGetNodePtr(void);
  int MemDebugGetNCountPtr(void);
#endif
};

template<class NodeType> class Tree {

private:

protected:

  // the vector containing the nodes
  std::vector<Pointer<NodeType> > nodes;
  
public:
  
  // constructor
  Tree();
  
  // destructor
  ~Tree();

  // Tree's methods
  
  // add a subnode
  void AddNode(const Pointer<NodeType> n);
  void AddCopyOfNode(const Pointer<NodeType> n);
  // delete a subnode
  bool DeleteNode(const Int i);
  // I would love to just return the iterators to this vector,
  // but GCC3.2 issues warnings against it, says it's deprecated.
  // I really don't see how, but still... Never mind, do it this
  // way. It's only ever used in Expression::ReorderNodes() anyway
  std::vector<Pointer<NodeType> >* GetNodeVectorPtr(void);
  // delete all subnodes
  void DeleteAllNodes(void);
  // get a subnode 
  Pointer<NodeType> GetNode(const Int i) const;
  // get a subnode 
  Pointer<NodeType>* GetNodePtr(const Int i);
  // get a copy of subnode 
  Pointer<NodeType> GetCopyOfNode(const Int i);
  // get the size of nodes
  Int GetSize(void) const;
  // compare two trees
  bool operator == (const Tree<NodeType>& t) const;
};

#endif

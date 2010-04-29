/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        tree.cxx                                               *
* Source:      GNU C++                                                *
* Purpose:     tree construction methods                              *
* History:     010517 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#define RCS7 "$Id: tree.cxx,v 1.13 2005/09/08 23:37:23 liberti Exp liberti $"

#include<iostream>
#include<vector>
#include<string>
#include "tree.h"
#include "exceptions.h"

////////////// MEMCHECK debug ////////////////////
#ifdef MEMDEBUG
vector<pair<void*,pair<int,int*> > > memcheckdebug;
vector<pair<void*,pair<int,int*> > >::iterator memcheckdebugit;
pair<void*, pair<int, int*> > memcheckdebugpair;
int memcheckdebugcounter;
#endif
///////////// END MEMCHECK debug ////////////////

// constructors  
template<class NodeType> Pointer<NodeType>::Pointer() {
  using namespace std;
  node = new NodeType;
  ncount = new Int(1);
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: c1 (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}

template<class NodeType> Pointer<NodeType>::Pointer(NodeType& n) {
  using namespace std;
  node = new NodeType(n);
  ncount = new Int(1);
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: c2 (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}

template<class NodeType> Pointer<NodeType>::Pointer(double v) {
  using namespace std;
  node = new NodeType(v);
  ncount = new Int(1);
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: c3 (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}

template<class NodeType> 
Pointer<NodeType>::Pointer(double c, int vi, std::string vn) {
  using namespace std;
  node = new NodeType(c, vi, vn);
  ncount = new Int(1);
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: c4 (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}
 
template<class NodeType> Pointer<NodeType>::Pointer(bool notinitialized) {
  using namespace std;
  if (!notinitialized) {
    node = new NodeType;
    ncount = new Int(1);
  } else {
    node = NULL;
    ncount = NULL;
  }
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: c5 (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}
  
// copy constructor
template<class NodeType> 
void Pointer<NodeType>::SetTo(const Pointer<NodeType>& t) {
  using namespace std;
  if (node != t.node) {
    Destroy();
    node = t.node;
    ncount = t.ncount;
  }
  (*ncount)++;
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: setto (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}    

template<class NodeType> Pointer<NodeType>::Pointer(const Pointer<NodeType>& t){
  using namespace std;
  node = NULL;
  ncount = NULL;
  SetTo(t);
}

// copy assignment: does NOT copy, it just references.
template<class NodeType> 
Pointer<NodeType>& Pointer<NodeType>::operator=(const Pointer<NodeType>& t) {
  using namespace std;
  SetTo(t);
  return *this;
}

// copy factory
// 1. this is a copy of pointer
template<class NodeType> 
void Pointer<NodeType>::SetToCopyOf(const Pointer<NodeType>& t) {
  using namespace std;
  if (node != t.node || node == t.node && ncount && *ncount == 1) {
    Destroy();
  } else if (node == t.node) {
    assert(ncount);
    *ncount--;
  }
  // destroys generality to use the user-defined constructor to force
  // the copy, but what the heck, can't work miracles.
  node = new NodeType(t, true); 
  ncount = new Int(1);
#ifdef MEMDEBUG
  memcheckdebugpair.first = (void*) node;
  memcheckdebugpair.second.first = memcheckdebugcounter;
  memcheckdebugpair.second.second = ncount;
  if (ncount && *ncount == 1)
    memcheckdebug.push_back(memcheckdebugpair);
  cerr << "MEMCHECK: setcpof (" << memcheckdebugcounter 
       << "): node = " << node << "; ncount = " 
       << ncount;
  memcheckdebugcounter++;
  if (ncount) {
    cerr << "; *ncount = " << *ncount << endl;
  } else {
    cerr << endl;
  }
#endif
}
// 2. returns a copy of this
template<class NodeType> Pointer<NodeType> Pointer<NodeType>::Copy(void) const{
  using namespace std;
  Pointer<NodeType> ret(true); // uninitialized
  ret.SetToCopyOf(*this);
  return ret;
}

// destructor
template<class NodeType> void Pointer<NodeType>::Destroy(void) {
  using namespace std;
#ifdef MEMDEBUG
  bool found = false;
  for(memcheckdebugit = memcheckdebug.begin(); 
      memcheckdebugit != memcheckdebug.end();
      memcheckdebugit++) {
    if (memcheckdebugit->first == (void*) node) {
      found = true;
      break;
    }
  }
  int lineno = -1;
  if (found) {
    lineno = memcheckdebugit->second.first;
    if  (ncount && *ncount == 1) {
      memcheckdebug.erase(memcheckdebugit);
    }
  }
  cerr << "MEMCHECK: deleting (" << lineno 
       << "): node = " << node << "; ncount = " 
       << ncount;
  if (ncount) {
    cerr << "; *ncount = " << *ncount - 1;
  }
#endif
  if (ncount) {
    if (--(*ncount) == 0) {
      if (node) {
	delete node; 
	node = NULL;
      }
      delete ncount;
      ncount = NULL;
#ifdef MEMDEBUG
      cerr << "; deleted";
#endif
    }
  }
#ifdef MEMDEBUG
  cerr << endl;
#endif    
}

template<class NodeType> Pointer<NodeType>::~Pointer() {
  using namespace std;
  Destroy();
}  

// overload of ->
template<class NodeType> NodeType* Pointer<NodeType>::operator->() const {
  using namespace std;
  return node;
}

template<class NodeType> NodeType Pointer<NodeType>::GetPointee(void) const {
  using namespace std;
  return *node;
}

// check for equality
template<class NodeType> 
bool Pointer<NodeType>::operator==(const Pointer<NodeType>& t) {
  using namespace std;
  if (node == t.node) {
    // fast check
    return true;
  } else {
    // use the NodeType::operator==
    if (*node == *(t.node)) {
      return true;
    }
  }
  return false;
}

#ifdef MEMDEBUG
template<class NodeType> NodeType* Pointer<NodeType>::MemDebugGetNodePtr(void){
  using namespace std;
  return node;
}
template<class NodeType> int Pointer<NodeType>::MemDebugGetNCountPtr(void) {
  using namespace std;
  return ncount;
}
#endif

// constructor
template<class NodeType> Tree<NodeType>::Tree() { }
  
// destructor
template<class NodeType> Tree<NodeType>::~Tree() { 
  /* // SIGSEGVs -- investigate or bear the memleak
     int sz = nodes.size();
     for(int i = 0; i < sz; i++) {
     nodes[i].Destroy();
     }
  */
  DeleteAllNodes();
}

// Tree's methods
  
// add a subnode
template<class NodeType> 
void Tree<NodeType>::AddNode(const Pointer<NodeType> n) {
  using namespace std;
  nodes.push_back(n);
}

template<class NodeType> 
void Tree<NodeType>::AddCopyOfNode(const Pointer<NodeType> n) {
  using namespace std;
  nodes.push_back(n.Copy());
}

// delete a subnode
template<class NodeType> bool Tree<NodeType>::DeleteNode(const Int i) {
  using namespace std;
  if (i >= (Int) nodes.size())
    return false;
  else {
    nodes.erase(nodes.begin() + i);
    return true;
  }
}  

// I would love to just return the iterators to this vector,
// but GCC3.2 issues warnings against it, says it's deprecated.
// I really don't see how, but still... Never mind, do it this
// way. It's only ever used in Expression::ReorderNodes() anyway
template<class NodeType> 
std::vector<Pointer<NodeType> >* Tree<NodeType>::GetNodeVectorPtr(void) {
  using namespace std;
  return &nodes;
}

// delete all subnodes
template<class NodeType> void Tree<NodeType>::DeleteAllNodes(void) {
  using namespace std;
  nodes.erase(nodes.begin(), nodes.end());    
}

// get a subnode 
template<class NodeType> 
Pointer<NodeType> Tree<NodeType>::GetNode(const Int i) const {
  using namespace std;
  //cout << "Tree::GetNode(): " << i << ", " << nodes.size() << endl;
  if (nodes.size() == 0) {
    throw ErrNotPermitted(30, "Tree", "GetNode", "nodes.size()==0",
			  "vector of nodes is empty", HELPURL);
  }
  if (i < (Int)nodes.size()) 
    return nodes[i];
  else
    throw ErrNotPermitted(32, "Tree", "GetNode", "i>=nodes.size()",
			  "not enough nodes in vector", HELPURL);
}

// get a subnode 
template<class NodeType> 
Pointer<NodeType>* Tree<NodeType>::GetNodePtr(const Int i) {
  using namespace std;
  //cout << "Tree::GetNode(): " << i << ", " << nodes.size() << endl;
  if (nodes.size() == 0) {
    throw ErrNotPermitted(33, "Tree", "GetNodePtr", "nodes.size()==0",
			  "vector of nodes is empty", HELPURL);
  }
  if (i < (Int)nodes.size()) 
    return &(nodes[i]);
  else
    throw ErrNotPermitted(34, "Tree", "GetNodePtr", "i>=nodes.size()",
			  "not enough nodes in vector", HELPURL);
}

// get a copy of subnode 
template<class NodeType> 
Pointer<NodeType> Tree<NodeType>::GetCopyOfNode(const Int i) {    
  using namespace std;
  Int j = i;
  if (nodes.size() == 0) {
    throw ErrNotPermitted(31, "Tree", "GetNode", "nodes.size()==0",
			  "vector of nodes is empty", HELPURL);
  }
  if (j >= (Int) nodes.size()) {
    throw ErrNotPermitted(35, "Tree", "GetCopyOfNode", "i>=nodes.size()",
			  "not enough nodes in vector", HELPURL);
  }
  Pointer<NodeType> ret;
  ret.SetToCopyOf(nodes[j]);
  return ret;
}

// get the size of nodes
template<class NodeType> Int Tree<NodeType>::GetSize(void) const {
  using namespace std;
  return nodes.size();
}

// compare two trees
template<class NodeType> 
bool Tree<NodeType>::operator==(const Tree<NodeType>& t) const {
  using namespace std;
  if (this == &t) {
    // fast check
    return true;
  } else {
    // recurse
    Int s = GetSize();
    if (s == t.GetSize()) {
      Int i;
      for(i = 0; i < s; i++) {
	if (!(GetNode(i) == t.GetNode(i)))
	  return false;
      }
      return true;
    } else
      return false;
  }
}




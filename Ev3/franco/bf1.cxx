/**********************************************************************
* Author:      Franco Raimondi                                        *
* Name:        bf.h                                                   *
* Source:      GNU C++                                                *
* Purpose:     Boolean expressions (base classes and functionality)   *
* History:     040609 work started                                    *
***********************************************************************/


using namespace std;

#include <string>
#if __GNUC__ == 3
#include <sstream>
#else
#include <strstream>
#endif
#include <cmath>
#include <algorithm>
#include <cassert>
#include "bf.h"

#define PEV3PI 3.1415926535897932385
#define NOTVARNAME "_var_not_found_"


///////////// classes ////////////////

// constructors
Operand::Operand() : 
  oplabel(CONST), constant(0), varindex(NOVARIABLE) { }


Operand::Operand(bool t) : 
  oplabel(CONST), constant(t), varindex(NOVARIABLE) { }


// create a variable leaf
Operand::Operand(const int t, string vn)
{ 
  // make it a variable
  oplabel = VAR;
  constant = 0;
  varindex = t;
  varname = vn;
}

// copy constructor
Operand::Operand(const Operand& t) : 
  oplabel(t.oplabel), constant(t.constant), 
  varindex(t.varindex), varname(t.varname) { }

// destructor
Operand::~Operand() { }

// Operand class methods:

string Operand::ToString(void) const {
#if __GNUC__ == 3
  stringstream outbuf;  
#else
  strstream outbuf;
#endif
  string vn;
  if (GetOpType() == CONST) {
    // constant
    outbuf << GetValue();
  } else if (GetOpType() == VAR) {
    // variable
    vn = GetVarName();
    outbuf << vn;
  } else {
    // operand, don't print anything
    ;
  }
  return outbuf.str();
}

// get operator type
int Operand::GetOpType(void) const { return oplabel; }
  
// get constant value - in CONSTs, multiply by coeff. and raise 
// to exponent, first
bool Operand::GetValue(void) const { 
  if (oplabel == CONST) {
    return constant;
  } else {
    return constant; 
  }
}

// get variable index
Int Operand::GetVarIndex(void) const { return varindex; }

// get variable name
string Operand::GetVarName(void) const { return varname; }

// set operator type
void Operand::SetOpType(const int t) { oplabel = t; }
  
// set constant value
void Operand::SetValue(const bool t) { oplabel = CONST; constant = t; }

// set variable index
void Operand::SetVarIndex(const Int t) { oplabel = VAR; varindex = t; }

// set variable name
void Operand::SetVarName(const string vn) { oplabel = VAR; varname = vn; }

// is operand a constant?
bool Operand::IsConstant(void) const {
  return (GetOpType() == CONST);
}

// is operand a variable?
bool Operand::IsVariable(void) const {
  return (GetOpType() == VAR);
}

// is operand a leaf node?
bool Operand::IsLeaf(void) const {
  return (IsConstant() || IsVariable());
}

// is operand a zero constant?
bool Operand::IsZero(void) const {
  if (GetOpType() == CONST) {
    bool c = GetValue();
    if (c) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
  

// is operand this == operand t?
bool Operand::operator == (const Operand& t) {
  if (this == &t) {
    // first fast check
    return true;
  } else {
    // not the same operand - check data fields
    return (GetOpType() == t.GetOpType() &&
	    GetValue() == t.GetValue() &&
	    GetVarIndex() == t.GetVarIndex());
  }
}      

// substitute a variable with a constant
void Operand::SubstituteVariableWithConstant(long int varindex, bool c) {
  if (GetOpType() == VAR && GetVarIndex() == varindex) {  
    SetOpType(CONST);
    SetVarIndex(NOVARIABLE);
    SetValue(c);    
  }
}


// create empty
BasicBoolForm::BasicBoolForm() { }

// create a constant leaf
BasicBoolForm::BasicBoolForm(const bool t) : Operand(t) { }

// create a variable leaf
BasicBoolForm::BasicBoolForm(const int t, const string vn) : 
  Operand(t, vn) { }
  
// user-defined copy constructor with two options
BasicBoolForm::BasicBoolForm(const BoolForm& t, const bool iscopy) :
  Operand(t.GetPointee()) {
  int i = 0;
  int s = t->GetSize();
  if (iscopy) {
    // create a _copy_ of t, subnode by subnode
    for ( ; i < s; i++)
      nodes.push_back(t->GetCopyOfNode(i));
  } else {
    // create a copy of the pointers in t
    for ( ; i < s; i++)
      nodes.push_back(t->GetNode(i));
  }
}

// copy constructor
BasicBoolForm::BasicBoolForm(const BasicBoolForm& t) : Operand(t) { 
  Int i = 0;
  Int s = t.GetSize();
  // necessary for constructs "(BasicBoolForm) e =" where e already 
  // existed prior to assignment  
  DeleteAllNodes(); 
  // create a copy of the subnode pointers in t
  for ( ; i < s; i++)
    nodes.push_back(t.GetNode(i));
}

// destructor
BasicBoolForm::~BasicBoolForm() { }

// BasicBoolForm class methods:

void BasicBoolForm::Debug (void) const {
  Int s = GetSize();
  cerr << "BasicBoolForm: Debug:\n";
  cerr << "\tthis   = " << this << endl;
  cerr << "\toptype = " << GetOpType() << endl;
  cerr << "\tnodes  = " << s << endl;
  Int i;
  for (i = 0; i < s; i++) {
    cerr << "\tnode " << i << ": " << GetNode(i)->GetOpType() << endl;
  }
}  

void BasicBoolForm::Zero(void) {
  DeleteAllNodes();
  SetValue(0);
  SetOpType(CONST);
}

void BasicBoolForm::One(void) {
  DeleteAllNodes();
  SetValue(1);
  SetOpType(CONST);
}

string BasicBoolForm::PrintTree(int blanks, int tabs) const {
#if __GNUC__ == 3
  stringstream outbuf;
#else
  strstream outbuf;
#endif 
  string b(blanks, ' ');
  if (IsLeaf()) {
    outbuf << b << Operand::ToString();
  } else {
    outbuf << b << "OP[" << GetOpType() << "](" << endl;
    for(int i = 0; i < GetSize(); i++) {
      outbuf << GetNode(i)->PrintTree(blanks + tabs, tabs);
      if (i < GetSize() - 1) {
        outbuf << ",";
      }
      outbuf << endl;
    }
    outbuf << b << ")";
  }
  return outbuf.str();
}

string BasicBoolForm::ToString(void) const {
#if __GNUC__ == 3
  stringstream outbuf;
#else
  strstream outbuf;
#endif  
  Int i, s;
  if (IsLeaf()) {
    // leaf node, use Operand::ToString()
    return Operand::ToString();
  } else {
    s = GetSize();
    if (s > 1) {
      string t;
      for (i = 0; i < s; i++) {
	t = GetNode(i)->ToString();
	outbuf << "(" << t << ")";
	if (i < s - 1) {
	  switch(GetOpType()) {
	  case AND:
	    outbuf << "*";
	    break;
	  case OR:
	    outbuf << "+";
	    break;
	  default:
	    outbuf << "UNKNOWNOP";
	    break;
	  }
	}
      }    
    } else {
      switch(GetOpType()) {
      case NOT:
	outbuf << "!";
	break;
      default:
	outbuf << "UNKNOWNOP";
      }
      if (s == 1) {
	string t(GetNode(0)->ToString());
	outbuf << "(" << t << ")";
      } else {
	// no arguments - error
	outbuf << "(NOARG)";
      }
    }
    return outbuf.str();
  }
}



void BasicBoolForm::VariableToConstant(long int varindex, bool c) {
  if (IsLeaf()) {
    // nonrecursive
    SubstituteVariableWithConstant(varindex, c);
  } else {
    // recursive
    int i;
    for (i = 0; i < GetSize(); i++) {
      GetNode(i)->VariableToConstant(varindex, c);
    }
  }
}



// OR operator:
BoolForm operator + (BoolForm a, BoolForm b) {
  BoolForm ret;

  if (a->IsLeaf() && a->GetOpType() == CONST &&
      b->IsLeaf() && b->GetOpType() == CONST) {
    // a, b are constants
    ret.SetToCopyOf(a);
    ret->SetValue(a->GetValue() || b->GetValue());
    return ret;
  } else {
    ret->SetOpType(OR);
    ret->AddCopyOfNode(a);
    ret->AddCopyOfNode(b);
    return ret;
  }
}

// AND operator:
BoolForm operator * (BoolForm a, BoolForm b) {
  BoolForm ret;

  if (a->IsLeaf() && a->GetOpType() == CONST &&
      b->IsLeaf() && b->GetOpType() == CONST) {
    // a, b are constants
    ret.SetToCopyOf(a);
    ret->SetValue(a->GetValue() && b->GetValue());
    return ret;
  } else {
    ret->SetOpType(AND);
    ret->AddCopyOfNode(a);
    ret->AddCopyOfNode(b);
    return ret;
  }
}

// NOT operator
BoolForm operator ! (BoolForm a) {
  BoolForm ret;
    
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    ret->SetValue(!(a->GetValue()));
  } else {
    ret->SetOpType(NOT);
    ret->AddCopyOfNode(a);
    
  }
  return ret;
}

// Push negations to atoms, using de Morgan
BoolForm PushNeg ( BoolForm a ) 
{
  BoolForm ret;
  int s;
  if ( a->IsLeaf() ) {
    ret.SetToCopyOf(a);
    return ret;
  } else {
    s = a->GetSize();
    if (s > 1) {
      BoolForm f1, f2;
      f1 = a->GetNode(0);
      f2 = a->GetNode(1);
      switch( a->GetOpType() ) {
      case AND:
	ret = PushNeg(f1) * PushNeg(f2);
	break;
      case OR:
	ret = PushNeg(f1) + PushNeg(f2);
	break;
      default:
	break;
      }
    } else {
      if (a->GetOpType() == NOT ) {
	// Difficult bit.
	BoolForm f1;
	f1 = a->GetNode(0);
	if ( f1->IsLeaf() ) {
	  ret.SetToCopyOf(a);
	  return ret;
	}
	else {
	  switch( f1->GetOpType() ) {
	  case AND:
	    { BoolForm f2,f3;
	    f2 = f1->GetNode(0);
	    f3 = f1->GetNode(1);	    
	    ret = PushNeg(!f2) + PushNeg(!f3);
	    return ret;
	    break;
	    }
	    
	  case OR:
	    { BoolForm f2,f3;
	    f2 = f1->GetNode(0);
	    f3 = f1->GetNode(1);
	    ret = PushNeg(!f2) * PushNeg(!f3);
	    return ret;
	    break; 
	    }
	  case NOT:
	    { BoolForm f2;
	    f2 = f1->GetNode(0);
	    ret = PushNeg(f2);
	    return ret;
	    break;
	    }
	    
	  default:
	    break;
	  }
	}
      }
    }
  }
  return ret;
}

  

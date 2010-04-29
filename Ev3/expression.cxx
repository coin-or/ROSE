/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        expression.cxx                                         *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (base classes and functionality)   *
* History:     010517 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#define RCS4 "$Id: expression.cxx,v 1.27 2006/07/30 05:36:39 liberti Exp liberti $"

#include <string>
#ifdef HAS_STRINGSTREAM
#include <sstream>
#else
#include <strstream>
#endif
#include <cmath>
#include <algorithm>
#include <cassert>
#include "expression.h"

using namespace std;

#define PEV3PI 3.1415926535897932385
#define NOTVARNAME "_var_not_found_"
#define VNAMEIDXCHAR "_"
//#define VNAMEIDXCHAR ""

///////////// BasicExpression class ////////////////

// create empty
BasicExpression::BasicExpression() { FastEvalAllocated = false; }

// create a constant leaf
BasicExpression::BasicExpression(const double t) : Operand(t) {
  FastEvalAllocated = false;
}

// create a constant (integer-valued) leaf
BasicExpression::BasicExpression(const Int t) : Operand(t) {
  FastEvalAllocated = false;
}

// create an (empty) operator or a variable leaf
BasicExpression::BasicExpression(const Int t, const bool isvar) :
  Operand(t, isvar) {
  FastEvalAllocated = false;
}

// create a variable leaf and set coefficient
BasicExpression::BasicExpression(const double c, const Int t, string vn) :
  Operand(c, t, vn) {
  FastEvalAllocated = false;
}

// user-defined copy constructor with two options
BasicExpression::BasicExpression(const Expression& t, const bool iscopy) :
  Operand(t.GetPointee()) {
  Int s = t->GetSize();
  if (iscopy) {
    // create a _copy_ of t, subnode by subnode
    for (int i = 0; i < s; i++) {
      nodes.push_back(t->GetCopyOfNode(i));
    }
  } else {
    // create a copy of the pointers in t
    for (int i = 0; i < s; i++) {
      nodes.push_back(t->GetNode(i));
    }
  }
  FastEvalAllocated = false;
}

// copy constructor
BasicExpression::BasicExpression(const BasicExpression& t) : Operand(t) {
  Int s = t.GetSize();
  // necessary for constructs "(BasicExpression) e =" where e already
  // existed prior to assignment
  // -- segvs (?)
  for(int i = 0; i < GetSize(); i++) {
    RecursiveDestroy(GetNodePtr(i));
  }
  DeleteAllNodes();
  // create a copy of the subnode pointers in t
  for (int i = 0; i < s; i++) {
    nodes.push_back(t.GetNode(i));
  }
  FastEvalAllocated = false;
}

// destructor
BasicExpression::~BasicExpression() {
  if (FastEvalAllocated) {
    DeleteFastEvalTree();
  }
}

// BasicExpression class methods:

int BasicExpression::GetNumberOfNodes(void) {
  int sz = GetSize();
  int ret = sz;
  if (!IsLeaf()) {
    for(int i = 0; i < sz; i++) {
      ret += GetNode(i)->GetNumberOfNodes();
    }
  }
  return ret;
}

int BasicExpression::GetNumberOfLevels(void) {
  int maxlevel = 1;
  GetNumberOfLevelsRecursive(maxlevel);
  return maxlevel;
}
void BasicExpression::GetNumberOfLevelsRecursive(int& maxlevel) {
  if (!IsLeaf()) {
    int thelevel = 0;
    int sz = GetSize();
    for(int i = 0; i < sz; i++) {
      int lvs = GetNode(i)->GetNumberOfLevels();
      if (lvs > thelevel) {
	thelevel = lvs;
      }
    }
    maxlevel += thelevel;
  }
}

# define DAGTOLERANCE 1e-8
int BasicExpression::GetNumberOfExplicitLevels(void) {
  int maxlevel = 1;
  GetNumberOfLevelsRecursive(maxlevel);
  return maxlevel;
}
void BasicExpression::GetNumberOfExplicitLevelsRecursive(int& maxlevel) {
  if (!IsLeaf()) {
    int thelevel = 0;
    int sz = GetSize();
    for(int i = 0; i < sz; i++) {
      int lvs = GetNode(i)->GetNumberOfLevels();
      if ( (fabs(GetCoeff() - 1) > DAGTOLERANCE &&
	    fabs(GetExponent() - 1) <= DAGTOLERANCE) ||
	   (fabs(GetCoeff() - 1) <= DAGTOLERANCE &&
	    fabs(GetExponent() - 1) > DAGTOLERANCE) ) {
	lvs += 2;
      } else if (fabs(GetCoeff() - 1) > DAGTOLERANCE &&
		 fabs(GetExponent() - 1) > DAGTOLERANCE) {
	lvs += 3;
      }
      if (lvs > thelevel) {
	thelevel = lvs;
      }
    }
    maxlevel += thelevel;
  }
}

int BasicExpression::GetNumberOfDistinctConstants(void) {
  map<double,int,DoubleLessThan> co;
  GetNumberOfDistinctConstantsRecursive(co);
  return co.size();
}
void BasicExpression::GetNumberOfDistinctConstantsRecursive
           (map<double,int,DoubleLessThan>& ConstOccurrence) {
  if (IsConstant()) {
    ConstOccurrence[GetValue()]++;
  } else {
    if (fabs(GetCoeff() - 1) > DAGTOLERANCE) {
      ConstOccurrence[GetCoeff()]++;
    }
    if (fabs(GetExponent() - 1) > DAGTOLERANCE) {
      ConstOccurrence[GetExponent()]++;
    }
    int sz = GetSize();
    for(int i = 0; i < sz; i++) {
      GetNode(i)->GetNumberOfDistinctConstantsRecursive(ConstOccurrence);
    }
  }
}

void BasicExpression::DeleteFastEvalTree(void) {
  if (FastEvalAllocated) {
    DeleteFastEvalTreeRecursive(fetree);
  }
  delete fetree;
  FastEvalAllocated = false;
}

void BasicExpression::DeleteFastEvalTreeRecursive(FastEvalTree* fet) {
  for(Int i = 0; i < fet->nodesize; i++) {
    DeleteFastEvalTreeRecursive(&fet->nodes[i]);
  }
  delete [] fet->nodes;
}

void BasicExpression::CreateFastEvalTree(void) {
  if (FastEvalAllocated) {
    DeleteFastEvalTree();
  }
  if (!FastEvalAllocated) {
    fetree = new FastEvalTree;
    CreateFastEvalTreeRecursive(fetree);
  }
}

void BasicExpression::CreateFastEvalTreeRecursive(FastEvalTree* fet) {
  Int s = GetSize();
  fet->optype = GetOpType();
  fet->coeff = GetCoeff();
  fet->depcoeff = depcoefficient;
  fet->exponent = GetExponent();
  fet->depexponent = depexponent;
  fet->value = GetValue();
  fet->depvalue = depconstant;
  fet->varindex = GetVarIndex();
  fet->nodesize = s;
  fet->nodes = new FastEvalTree [fet->nodesize];
  for(Int i = 0; i < s; i++) {
    GetNode(i)->CreateFastEvalTreeRecursive(&fet->nodes[i]);
  }
}

void BasicExpression::Debug (void) const {
  Int s = GetSize();
  cerr << "BasicExpression: Debug:\n";
  cerr << "\tthis   = " << this << endl;
  cerr << "\toptype = " << GetOpType() << endl;
  cerr << "\tnodes  = " << s << endl;
  Int i;
  for (i = 0; i < s; i++) {
    cerr << "\tnode " << i << ": " << GetNode(i)->GetOpType() << endl;
  }
}

void BasicExpression::Zero(void) {
  // -- segvs (?)
  for(int i = 0; i < GetSize(); i++) {
    RecursiveDestroy(GetNodePtr(i));
  }
  DeleteAllNodes();
  SetCoeff(1.0);
  SetExponent(1.0);
  SetValue(0.0);
  SetOpType(CONST);
}

void BasicExpression::One(void) {
  // -- segvs (?)
  for(int i = 0; i < GetSize(); i++) {
    RecursiveDestroy(GetNodePtr(i));
  }
  DeleteAllNodes();
  SetCoeff(1.0);
  SetExponent(1.0);
  SetValue(1.0);
  SetOpType(CONST);
}

string BasicExpression::PrintTree(int blanks, int tabs) const {
#ifdef HAS_STRINGSTREAM
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

string BasicExpression::ToString(void) const {
#ifdef HAS_STRINGSTREAM
  stringstream outbuf;
#else
  strstream outbuf;
#endif
  Int i, s;
  if (IsLeaf()) {
    // leaf node, use Operand::ToString()
    return Operand::ToString();
  } else {
    // operator node
    if (GetCoeff() != 1) {
      outbuf << "(" << GetCoeff() << "*(";
    }
    s = GetSize();
    if (s > 1) {
      string t;
      for (i = 0; i < s; i++) {
	t = GetNode(i)->ToString();
	outbuf << "(" << t << ")";
	if (i < s - 1) {
	  switch(GetOpType()) {
	  case SUM:
	    outbuf << "+";
	    break;
	  case DIFFERENCE:
	    outbuf << "-";
	    break;
	  case PRODUCT:
	    outbuf << "*";
	    break;
	  case FRACTION:
	    outbuf << "/";
	    break;
	  case POWER:
	    outbuf << "^";
	    break;
	  default:
	    outbuf << "UNKNOWNOP";
	    break;
	  }
	}
      }
    } else {
      switch(GetOpType()) {
      case PLUS:
	break;
      case MINUS:
	outbuf << "-";
	break;
      case LOG:
	outbuf << "log";
	break;
      case EXP:
	outbuf << "exp";
	break;
      case SIN:
	outbuf << "sin";
	break;
      case COS:
	outbuf << "cos";
	break;
      case TAN:
	outbuf << "tan";
	break;
      case COT:
	outbuf << "cot";
	break;
      case SINH:
	outbuf << "sinh";
	break;
      case COSH:
	outbuf << "cosh";
	break;
      case TANH:
	outbuf << "tanh";
	break;
      case COTH:
	outbuf << "coth";
	break;
      case SQRT:
	outbuf << "sqrt";
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
    if (GetCoeff() != 1) {
      outbuf << "))";
    }
    return outbuf.str();
  }
}

// is expression this == expression t?
bool BasicExpression::IsEqualTo(const Expression& t) const {
  //Just in case constant value not consolidated in coeff
  if (GetOpType() == CONST || t->GetOpType() == CONST) {
    if (GetValue() == t->GetValue()) {
      return true;
    }
    else {
      return false;
    }
  }
  
  if (IsEqualToNoCoeff(t)) {
    if (GetCoeff() == t->GetCoeff())
      return true;
    else
      return false;
  } else
    return false;
}

bool BasicExpression::operator == (const BasicExpression& t) const {
  // fast checks
  if (IsLeaf() && t.IsLeaf()) {
    if (GetOpType() == t.GetOpType()) {
      if (GetOpType() == CONST) {
	// if both are constants, they're always equal up to coefficient
	return true;
      } else if (GetOpType() == VAR && GetVarIndex() == t.GetVarIndex() &&
		 GetExponent() == t.GetExponent())
	return true;
      else
	return false;
    } else
      return false;
  } else if ((!IsLeaf()) && (!t.IsLeaf())) {
    // both BasicExpressions are not leaves, recurse using PRECISE
    // (i.e. not up to coefficient) form
    if (GetSize() != t.GetSize())
      return false;
    if (GetOpType() != t.GetOpType())
      return false;
    Int i;
    bool ret = true;
    for (i = 0; i < GetSize(); i++) {
      if (!(GetNode(i)->IsEqualTo(t.GetNode(i)))) {
	ret = false;
	break;
      }
    }
    if (ret) {
      return true;
    }
    return false;
  } else
    return false;
}

bool BasicExpression::IsEqualTo(const double t) const {
  return (IsLeaf() && GetOpType() == CONST && GetValue() == t);
}

int BasicExpression::NumberOfVariables() const {
  int maxvi = 0;
  return NumberOfVariables(maxvi);
}

int BasicExpression::NumberOfVariables(int& maxvi) const {
  int newvi = 0;
  if (IsVariable()) {
    newvi = GetVarIndex();
    if (newvi > maxvi) {
      maxvi = newvi;
    }
    return maxvi;
  } else if (!IsLeaf()) {
    for(int i = 0; i < GetSize(); i++) {
      newvi = GetNode(i)->NumberOfVariables(maxvi);
      if (newvi > maxvi) {
	maxvi = newvi;
      }
    }
    return maxvi;
  } else {
    return 0;
  }
}


// taken and adapted from operator ==
bool BasicExpression::IsEqualToNoCoeff(const Expression& t) const {
  // fast checks
  if (IsLeaf() && t->IsLeaf()) {
    if (GetOpType() == t->GetOpType()) {
      if (GetOpType() == CONST) {
	// if both are constants, they're always equal up to coefficient
	return true;
      } else if (GetOpType() == VAR && GetVarIndex() == t->GetVarIndex() &&
		 GetExponent() == t->GetExponent())
	return true;
      else
	return false;
    } else
      return false;
  } else if ((!IsLeaf()) && (!t->IsLeaf())) {
    // both BasicExpressions are not leaves, recurse using PRECISE
    // (i.e. not up to coefficient) form
    if (GetSize() != t->GetSize() || GetOpType() != t->GetOpType() ||
	GetExponent() != t->GetExponent()) {
      return false;
    }
    Int i;
    bool ret = true;
    for (i = 0; i < GetSize(); i++) {
      if (!(GetNode(i)->IsEqualTo(t->GetNode(i)))) {
	ret = false;
	break;
      }
    }
    if (ret) {
      return true;
    }
    return false;
  } else
    return false;
}

bool BasicExpression::IsEqualBySchema(const Expression& t) const {
  if (IsLeaf() && t->IsLeaf()) {
    if (GetOpType() == t->GetOpType() &&
	IsLinear() == t->IsLinear()) {
      return true;
    } else {
      return false;
    }
  } else if ((!IsLeaf()) && (!t->IsLeaf())) {
    // both BasicExpressions are not leaves, recurse
    if (GetSize() != t->GetSize()) {
      return false;
    }
    if (GetOpType() != t->GetOpType()) {
      return false;
    }
    Int i;
    bool ret = true;
    for (i = 0; i < GetSize(); i++) {
      if (!(GetNode(i)->IsEqualBySchema(t->GetNode(i)))) {
	ret = false;
	break;
      }
    }
    if (ret) {
      return true;
    }
    return false;
  } else {
    return false;
  }
}

bool BasicExpression::IsEqualByOperator(const int theoplabel) const {
  if (GetOpType() == theoplabel) {
    return true;
  } else {
    return false;
  }
}

double BasicExpression::Eval(double* varvalues, Int vsize) const {
  map<Int,Int> varmap;
  for(Int i = 1; i <= vsize; i++) {
    varmap[i] = i;
  }
  return Eval(varvalues, varmap, vsize);
}

double BasicExpression::Eval(double* varvalues, map<int,int>& varmap,
			     Int vsize) const {
  // evaluate expression
  if (GetOpType() == CONST) {
    return GetValue();
  } else if (GetOpType() == VAR) {
    if (GetVarIndex() > vsize) {
      throw ErrNotPermitted(0, "BasicExpression", "Eval", "varindex > vsize",
			    "varindex should be <= vsize", HELPURL);
    } else if (GetVarIndex() <= 0) {
      throw ErrNotPermitted(1, "BasicExpression", "Eval", "varindex <= 0",
                            "varindex should be > 0", HELPURL);
    }
    double ret = GetCoeff() * pow(varvalues[varmap[varindex] - 1],
				  GetExponent());
    return ret;
  } else {
    Int i;
    double ret = 0;
    double tmp = 0;
    int op = GetOpType();
    if (GetSize() == 0) {
      throw ErrNotPermitted(5, "BasicExpression", "Eval", "GetSize()==0",
			    "non-leaf expressions without subnodes",
			    HELPURL);
    }
    switch(op) {
    case SUM:
      for(i = 0; i < GetSize(); i++) {
	ret += GetNode(i)->Eval(varvalues, varmap, vsize);
      }
      ret *= GetCoeff();
      break;
    case DIFFERENCE:
      for(i = 0; i < GetSize(); i++) {
	ret -= GetNode(i)->Eval(varvalues, varmap, vsize);
      }
      ret *= GetCoeff();
      break;
    case PRODUCT:
      ret = 1;
      for(i = 0; i < GetSize(); i++) {
	ret *= GetNode(i)->Eval(varvalues, varmap, vsize);
      }
      ret *= GetCoeff();
      break;
    case FRACTION:
      if (GetSize() != 2) {
	cerr << "BasicExpression::Eval: in [GetSize()!=2]: "
	     << "fractions should have just 2 operands, but going ahead "
	     << "anyway...\n";
	//throw ErrNotPermitted(1, "BasicExpression", "Eval", "GetSize()!=2",
	//"fractions should have just 2 operands",
	//HELPURL);
      }
      for(i = 0; i < GetSize(); i++) {
	tmp = GetNode(i)->Eval(varvalues, varmap, vsize);
	if (i > 0 && tmp == 0) {
	  cerr << "BasicExpression::Eval: division by zero not allowed (node "
	       << i << " evaluates to 0), setting to a large value" << endl;
	  //throw ErrDivideByZero(2, "BasicExpression", "Eval", "divisor==0",
	  //"division by zero not allowed", HELPURL);
	  tmp = Ev3NearZero();
	}
	if (i == 0)
	  ret = tmp;
	else
	  ret /= tmp;
      }
      ret *= GetCoeff();
      break;
    case POWER:
      if (GetSize() != 2) {
	cerr << "BasicExpression::Eval: in [GetSize()!=2]: "
	     << "powers should have just 2 operands, but going ahead "
	     << "anyway...\n";
	//throw ErrNotPermitted(3, "BasicExpression", "Eval", "GetSize()!=2",
	//"powers should have just 2 operands",
	//HELPURL);
      }
      ret = GetNode(0)->Eval(varvalues, varmap, vsize);
      for(i = 1; i < GetSize(); i++) {
	ret = pow(ret, GetNode(i)->Eval(varvalues, varmap, vsize));
      }
      ret *= GetCoeff();
      break;
    // for all unary functions - should check that GetSize() == 1
    case PLUS:
      ret = GetCoeff() * GetNode(0)->Eval(varvalues, varmap, vsize);
      break;
    case MINUS:
      ret = -GetCoeff() * GetNode(0)->Eval(varvalues, varmap, vsize);
      break;
    case LOG:
      tmp = GetNode(0)->Eval(varvalues, varmap, vsize);
      if (tmp < 0) {
	cerr << "BasicExpression::Eval: log of negative not allowed ("
	     << "argument evaluates to < 0), taking absolute value" << endl;
	//throw ErrNotPermitted(6, "BasicExpression", "Eval", "log arg < 0",
	//"log of negative not allowed", HELPURL);
	tmp = -tmp;
      } else if (tmp == 0) {
	cerr << "BasicExpression::Eval: log of zero not allowed ("
	     << "argument evaluates to 0), setting to large negative value"
	     << endl;
	//throw ErrNotPermitted(7, "BasicExpression", "Eval", "log arg == 0",
	//"log of zero not allowed", HELPURL);
	tmp = Ev3NearZero();
      }
      ret = GetCoeff() * log(tmp);
      break;
    case EXP:
      ret = GetCoeff() * exp(GetNode(0)->Eval(varvalues, varmap, vsize));
      break;
    case SIN:
      ret = GetCoeff() * sin(GetNode(0)->Eval(varvalues, varmap, vsize));
      break;
    case COS:
      ret = GetCoeff() * cos(GetNode(0)->Eval(varvalues, varmap, vsize));
      break;
    case TAN:
      ret = GetCoeff() * tan(GetNode(0)->Eval(varvalues, varmap, vsize));
      // should check that cos() is nonzero
      break;
    case SINH:
      ret = GetCoeff() * sinh(GetNode(0)->Eval(varvalues, varmap, vsize));
      break;
    case COSH:
      ret = GetCoeff() * cosh(GetNode(0)->Eval(varvalues, varmap, vsize));
      break;
    case TANH:
      ret = GetCoeff() * tanh(GetNode(0)->Eval(varvalues, varmap, vsize));
      // should check that cosh is nonzero
      break;
    case COTH:
      ret = GetCoeff() / tanh(GetNode(0)->Eval(varvalues, varmap, vsize));
      // should check that tanh is nonzero
      break;
    case SQRT:
      tmp = GetNode(0)->Eval(varvalues, varmap, vsize);
      if (tmp < 0) {
	cerr << "BasicExpression::Eval: sqrt of negative not allowed, "
	     << "taking absolute value" << endl;
	//throw ErrNotPermitted(9, "BasicExpression", "Eval", "sqrt arg < 0",
	//"sqrt of negative not allowed", HELPURL);
	tmp = -tmp;
      }
      ret = GetCoeff() * sqrt(tmp);
      break;
    }
    return ret;
  }
}

double BasicExpression::FastEval(double* varvalues, Int vsize) {
  if (!FastEvalAllocated) {
    CreateFastEvalTree();
  }
  return FastEvalRecursive(fetree, varvalues, NULL, vsize);
}

double BasicExpression::FastEval(double* varvalues, map<int,int>& varmap,
				 Int vsize) {
  if (!FastEvalAllocated) {
    CreateFastEvalTree();
  }
  return FastEvalRecursive(fetree, varvalues, &varmap, vsize);
}

double BasicExpression::FastEval(FastEvalTree* fet,
				 double* varvalues, Int vsize) {
  if (!FastEvalAllocated) {
    CreateFastEvalTree();
  }
  return FastEvalRecursive(fet, varvalues, NULL, vsize);
}

double BasicExpression::FastEval(FastEvalTree* fet, double* varvalues,
				 map<int,int>& varmap, Int vsize) {
  if (!FastEvalAllocated) {
    CreateFastEvalTree();
  }
  return FastEvalRecursive(fet, varvalues, &varmap, vsize);
}

double BasicExpression::FastEvalRecursive(FastEvalTree* fet,
					  double* varvalues,
					  map<int,int>* varmap,
					  Int vsize) {
  double thecoeff = fet->depcoeff ? *(fet->depcoeff) : fet->coeff;
  double theexpon = fet->depexponent ? *(fet->depexponent) : fet->exponent;
  double thevalue = fet->depvalue ? *(fet->depvalue) : fet->value;
  // evaluate expression
  if (fet->optype == CONST) {
    return thevalue;
  } else if (fet->optype == VAR) {
    if (fet->varindex > vsize) {
      throw ErrNotPermitted(0, "BasicExpression", "FastEvalRecursive",
			    "fet->varindex > vsize",
			    "fet->varindex should be <= vsize", HELPURL);
    } else if (fet->varindex <= 0) {
      throw ErrNotPermitted(1, "BasicExpression", "FastEvalRecursive",
			    "fet->varindex <= 0",
			    "fet->varindex should be > 0",
			    HELPURL);
    }
//
   //cout << " FET: " << fet->varindex << endl;
   //cout << " FET: " << varvalues[fet->varindex - 1] << endl;
    double ret = 0;
    if (varmap) {
      ret = varvalues[(*varmap)[fet->varindex] - 1];
    } else {
      ret = varvalues[fet->varindex - 1];
    }
    ret = thecoeff * pow(ret, theexpon);
    return ret;
  } else {
    Int i;
    double ret = 0;
    double tmp = 0;
    int op = fet->optype;
    if (fet->nodesize == 0) {
      throw ErrNotPermitted(5, "BasicExpression", "FastEvalRecursive",
			    "GetSize()==0",
			    "non-leaf expressions without subnodes",
			    HELPURL);
    }
    switch(op) {
    case SUM:
      for(i = 0; i < fet->nodesize; i++) {
	ret += FastEvalRecursive(&fet->nodes[i], varvalues, varmap, vsize);
      }
      ret *= thecoeff;
      break;
    case DIFFERENCE:
      for(i = 0; i < fet->nodesize; i++) {
	ret -= FastEvalRecursive(&fet->nodes[i], varvalues, varmap, vsize);
      }
      ret *= thecoeff;
      break;
    case PRODUCT:
      ret = 1;
      for(i = 0; i < fet->nodesize; i++) {
	ret *= FastEvalRecursive(&fet->nodes[i], varvalues, varmap, vsize);
      }
      ret *= thecoeff;
      break;
    case FRACTION:
      if (fet->nodesize != 2) {
	cerr << "BasicExpression::FastEvalRecursive: in [GetSize()!=2]: "
	     << "fractions should have just 2 operands, but going ahead "
	     << "anyway...\n";
	//throw ErrNotPermitted(1, "BasicExpression", "Eval", "GetSize()!=2",
	//"fractions should have just 2 operands",
	//HELPURL);
      }
      for(i = 0; i < fet->nodesize; i++) {
	tmp = FastEvalRecursive(&fet->nodes[i], varvalues, varmap, vsize);
	if (i > 0 && tmp == 0) {
	  cerr << "BasicExpression::FastEvalRecursive: "
	       << "division by zero not allowed (node "
	       << i << " evaluates to 0), setting to a large value" << endl;
	  //throw ErrDivideByZero(2, "BasicExpression", "Eval", "divisor==0",
	  //"division by zero not allowed", HELPURL);
	  tmp = Ev3NearZero();
	}
	if (i == 0)
	  ret = tmp;
	else
	  ret /= tmp;
      }
      ret *= thecoeff;
      break;
    case POWER:
      if (fet->nodesize != 2) {
	cerr << "BasicExpression::FastEvalRecursive: in [GetSize()!=2]: "
	     << "powers should have just 2 operands, but going ahead "
	     << "anyway...\n";
	//throw ErrNotPermitted(3, "BasicExpression", "Eval", "GetSize()!=2",
	//"powers should have just 2 operands",
	//HELPURL);
      }
      ret = FastEvalRecursive(&fet->nodes[0], varvalues, varmap, vsize);
      for(i = 1; i < fet->nodesize; i++) {
	ret = pow(ret, FastEvalRecursive(&fet->nodes[i], varvalues, varmap,
					 vsize));
      }
      ret *= thecoeff;
      break;
    // for all unary functions - should check that GetSize() == 1
    case PLUS:
      ret = thecoeff * FastEvalRecursive(&fet->nodes[0],
					   varvalues, varmap, vsize);
      break;
    case MINUS:
      ret = -thecoeff * FastEvalRecursive(&fet->nodes[0],
					    varvalues, varmap, vsize);
      break;
    case LOG:
      tmp = FastEvalRecursive(&fet->nodes[0], varvalues, varmap, vsize);
      if (tmp < 0) {
	cerr << "BasicExpression::FastEvalRecursive: log of negative not allowed ("
	     << "argument evaluates to < 0), taking absolute value" << endl;
	//throw ErrNotPermitted(6, "BasicExpression", "Eval", "log arg < 0",
	//"log of negative not allowed", HELPURL);
	tmp = -tmp;
      } else if (tmp == 0) {
	cerr << "BasicExpression::FastEvalRecursive: log of zero not allowed ("
	     << "argument evaluates to 0), setting to large negative value"
	     << endl;
	//throw ErrNotPermitted(7, "BasicExpression", "Eval", "log arg == 0",
	//"log of zero not allowed", HELPURL);
	tmp = Ev3NearZero();
      }
      ret = thecoeff * log(tmp);
      break;
    case EXP:
      ret = thecoeff * exp(FastEvalRecursive(&fet->nodes[0],
					       varvalues, varmap, vsize));
      break;
    case SIN:
      ret = thecoeff * sin(FastEvalRecursive(&fet->nodes[0],
					       varvalues, varmap, vsize));
      break;
    case COS:
      ret = thecoeff * cos(FastEvalRecursive(&fet->nodes[0],
					       varvalues, varmap, vsize));
      break;
    case TAN:
      ret = thecoeff * tan(FastEvalRecursive(&fet->nodes[0],
					       varvalues, varmap, vsize));
      // should check that cos() is nonzero
      break;
    case SINH:
      ret = thecoeff * sinh(FastEvalRecursive(&fet->nodes[0],
						varvalues, varmap, vsize));
      break;
    case COSH:
      ret = thecoeff * cosh(FastEvalRecursive(&fet->nodes[0],
						varvalues, varmap, vsize));
      break;
    case TANH:
      ret = thecoeff * tanh(FastEvalRecursive(&fet->nodes[0],
						varvalues, varmap, vsize));
      // should check that cosh is nonzero
      break;
    case COTH:
      ret = thecoeff / tanh(FastEvalRecursive(&fet->nodes[0],
					      varvalues, varmap, vsize));
      // should check that tanh is nonzero
      break;
    case SQRT:
      tmp = FastEvalRecursive(&fet->nodes[0], varvalues, varmap, vsize);
      if (tmp < 0) {
	cerr << "BasicExpression::FastEvalRecursive: "
	     << "sqrt of negative not allowed, "
	     << "taking absolute value" << endl;
	//throw ErrNotPermitted(9, "BasicExpression", "Eval", "sqrt arg < 0",
	//"sqrt of negative not allowed", HELPURL);
	tmp = -tmp;
      }
      ret = thecoeff * sqrt(tmp);
      break;
    }
    return ret;
  }
}

FastEvalTree* BasicExpression::GetFastEvalTree(void) {
  if (!FastEvalAllocated) {
    CreateFastEvalTree();
  }
  return fetree;
}

bool BasicExpression::DependsOnVariable(Int vi) const {
  if (IsLeaf()) {
    if (GetOpType() == VAR) {
      if (GetVarIndex() == vi)
	return true;
      else
	return false;
    } else
      return false;
  } else {
    Int i;
    bool ret = false;
    for (i = 0; i < GetSize(); i++) {
      ret = GetNode(i)->DependsOnVariable(vi);
      if (ret)
	return true;
    }
    return false;
  }
}

int BasicExpression::DependsLinearlyOnVariable(Int vi) const {
  if (IsVariable()) {
    if (GetVarIndex() == vi) {
      if (GetExponent() == 1) {
	return 1; // depends linearly
      } else {
	return 0; // depends nonlinearly
      }
    } else {
      return 2; // doesn't depend on vi at all
    }
  } else {
    int i;
    int d;
    bool dependsatall = false;
    // if node is linear:
    if (GetOpType() == SUM || GetOpType() == DIFFERENCE ||
	GetOpType() == PLUS || GetOpType() == MINUS) {
      for(i = 0; i < GetSize(); i++) {
	d = GetNode(i)->DependsLinearlyOnVariable(vi);
	if (d == 0) {
	  // depends nonlinearly, return 0
	  return 0;
	}
	if (d == 1) {
	  // depends linearly, record
	  dependsatall = true;
	}
      }
      if (dependsatall) {
	return 1;
      } else {
	return 2;
      }
    } else {
      if (DependsOnVariable(vi)) {
	return 0;  // depends nonlinearly
      } else {
	return 2;  // doesn't depend on vi at all
      }
    }
  }
}

void BasicExpression::ConsolidateProductCoeffs(void) {
  if (GetOpType() == PRODUCT) {
    Int i;
    double tc = GetCoeff();
    for (i = 0; i < GetSize(); i++) {
      tc *= GetNode(i)->GetCoeff();
      GetNode(i)->SetCoeff(1);
    }
    if (fabs(tc) < Ev3NearZero()) {
      Zero();
    } else {
      SetCoeff(tc);
    }
  }
}

void BasicExpression::DistributeCoeffOverSum(void) {
  if (GetOpType() == SUM) {
    double tc = GetCoeff();
    if (tc != 1) {
      SetCoeff(1);
      Int i;
      for(i = 0; i < GetSize(); i++) {
	GetNode(i)->SetCoeff(tc * GetNode(i)->GetCoeff());
	GetNode(i)->DistributeCoeffOverSum();
      }
    }
  }
}

void BasicExpression::DistributeCoeffOverProduct(void) {
  if (GetOpType() == PRODUCT) {
    double tc = GetCoeff();
    if (tc != 1 && GetSize() > 0) {
      SetCoeff(1);
      GetNode(0)->SetCoeff(tc * GetNode(0)->GetCoeff());
    }
  }
}

// enforce constant dependencies (added for MORON - see ../PROGNOTES)
// this only acts on the proper leaf nodes
void BasicExpression::EnforceDependency(void) {
  if (IsLeaf()) {
    // nonrecursive
    EnforceDependencyOnOperand();
  } else {
    // recursive
    int i;
    for (i = 0; i < GetSize(); i++) {
      GetNode(i)->EnforceDependency();
    }
  }
}

// substitute a variable with a constant
void BasicExpression::VariableToConstant(int varindex, double c) {
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

// replace variable indexed v1 with variable indexed v2
void BasicExpression::ReplaceVariable(int v1, int v2, string vn) {
  if (DependsOnVariable(v1)) {
    if (IsVariable() && GetVarIndex() == v1) {
      SetVarIndex(v2);
      SetVarName(vn);
    } else {
      int i;
      for(i = 0; i < GetSize(); i++) {
	GetNode(i)->ReplaceVariable(v1, v2, vn);
      }
    }
  }
}
void BasicExpression::ReplaceVariable(int v1, int v2, string vn,
				      double c2) {
  if (DependsOnVariable(v1)) {
    if (IsVariable() && GetVarIndex() == v1) {
      SetVarIndex(v2);
      SetVarName(vn);
      SetCoeff(GetCoeff() * c2);
    } else {
      int i;
      for(i = 0; i < GetSize(); i++) {
	GetNode(i)->ReplaceVariable(v1, v2, vn, c2);
      }
    }
  }
}

// replace with a variable the deepest node conforming to schema and
// return replaced term or zero expression if no replacement occurs
// vn is "w", var names are cat(vn, vi)
Expression BasicExpression::ReplaceBySchema(int vi, string vn,
					    Expression schema) {
  Expression ret(0.0);
  ret = ReplaceBySchemaRecursive(vi, vn, schema);
  if (ret->IsZero()) {
    // no subnodes with schema found, wor on this one
    if (IsEqualBySchema(schema)) {
      // this node is according to schema
      // save (recursively) this into ret
      ret.SetToCopyOf(*this);
      // replace with w_vi
      // -- segvs (?)
      for(int i = 0; i < GetSize(); i++) {
	RecursiveDestroy(GetNodePtr(i));
      }
      DeleteAllNodes();
      SetOpType(VAR);
      SetVarIndex(vi);
      std::stringstream ss;
      ss << vn << vi;
      SetVarName(ss.str());
      SetCoeff(1.0);
      SetExponent(1.0);
    }
  }
  return ret;
}

// recursive version - works on subnodes, not on current node
// vn is "w", var names are cat(vn, vi)
Expression BasicExpression::ReplaceBySchemaRecursive(int vi, string vn,
						     Expression schema) {
  bool done = false;
  Expression ret(0.0);
  for(int i = 0; i < GetSize(); i++) {
    if (!GetNode(i)->IsLeaf()) {
      ret = GetNode(i)->ReplaceBySchemaRecursive(vi, vn, schema);
      if (!ret->IsZero()) {
	done = true;
	break;
      }
    }
    if (!done) {
      if (GetNode(i)->IsEqualBySchema(schema)) {
	ret = GetNode(i);
	std::stringstream ss;
	ss << vn << vi;
	Expression w(1, vi, ss.str());
	GetNodePtr(i)->SetTo(w);
	done = true;
	break;
      }
    }
  }
  return ret;
}

// replace with a variable the deepest node with given operator label
// return replaced term or zero expression if no replacement occurs
// vn is "w", var names are cat(vn, vi)
Expression BasicExpression::ReplaceByOperator(int vi, string vn,
					      int theoplabel) {
  Expression ret(0.0);
  ret = ReplaceByOperatorRecursive(vi, vn, theoplabel);
  if (ret->IsZero()) {
    // no subnodes with schema found, wor on this one
    if (IsEqualByOperator(theoplabel)) {
      // this node is according to schema
      // save (recursively) this into ret
      ret.SetToCopyOf(*this);
      // replace with w_vi
      // -- segvs
      for(int i = 0; i < GetSize(); i++) {
	RecursiveDestroy(GetNodePtr(i));
      }
      DeleteAllNodes();
      SetOpType(VAR);
      SetVarIndex(vi);
      std::stringstream ss;
      ss << vn << vi;
      SetVarName(ss.str());
      SetCoeff(1.0);
      SetExponent(1.0);
    }
  }
  return ret;
}

// recursive version - works on subnodes, not on current node
// vn is "w", var names are cat(vn, vi)
Expression BasicExpression::ReplaceByOperatorRecursive(int vi, string vn,
						       int theoplabel) {
  bool done = false;
  Expression ret(0.0);
  for(int i = 0; i < GetSize(); i++) {
    if (!GetNode(i)->IsLeaf()) {
      ret = GetNode(i)->ReplaceByOperatorRecursive(vi, vn, theoplabel);
      if (!ret->IsZero()) {
	done = true;
	break;
      }
    }
    if (!done) {
      if (GetNode(i)->IsEqualByOperator(theoplabel)) {
	ret = GetNode(i);
	std::stringstream ss;
	ss << vn << vi;
	Expression w(1, vi, ss.str());
	GetNodePtr(i)->SetTo(w);
	done = true;
	break;
      }
    }
  }

  return ret;
}

void BasicExpression::ReplaceWithExpression(Expression replace) {
  /* // -- segvs (?)
  for(int i = 0; i < GetSize(); i++) {
    RecursiveDestroy(GetNodePtr(i));
  }
  */
  DeleteAllNodes();
  if (replace->GetOpType() == VAR) {
    SetVarIndex(replace->GetVarIndex());
    SetVarName(replace->GetVarName());
    SetCoeff(replace->GetCoeff());
    SetExponent(replace->GetExponent());
  } else if (replace->GetOpType() == CONST) {
    SetValue(replace->GetValue());
  } else {
    SetCoeff(replace->GetCoeff());
    SetExponent(replace->GetExponent());
    SetOpType(replace->GetOpType());
  }
  for (int i = 0; i < replace->GetSize(); i++) {
    nodes.push_back(replace->GetNode(i));
  }
}


int BasicExpression::ReplaceSubexpression(Expression needle,
					  Expression replace) {
  int ret = 0;
  if (!IsLeaf()) {
    // recurse
    for(int i = 0; i < GetSize(); i++) {
      ret += GetNode(i)->ReplaceSubexpression(needle, replace);
    }
  }
  // act on this node
  if (IsEqualTo(needle)) {
    ret++;
    ReplaceWithExpression(replace);
  }
  return ret;
}

// smith's standard form
// defvarname is "w", var names are cat(defvarname, vi)
void BasicExpression::SmithStandardForm(int defvarindex,
					int addedvarstartindex,
					string defvarname,
					vector<int>& oplabels,
					vector<Expression>& schemata,
					vector<Expression>& defcons) {
  SmithStandardFormRecursive(defvarindex, addedvarstartindex, defvarname,
			     oplabels, schemata, defcons);
}

int BasicExpression::SmithFindInDefcons(int addedvarstartindex,
					std::vector<Expression>& defcons) {
  int ret = -1;
  int dcounter = 0;
  for(vector<Expression>::iterator dcit = defcons.begin();
      dcit != defcons.end(); dcit++) {
    if (IsEqualTo(*dcit)) {
      ret = dcounter;
      break;
    }
    dcounter++;
  }
  if (ret >= 0) {
    ret += addedvarstartindex;
  }
  return ret;
}

// defvarname is "w", var names are cat(defvarname, vi)
void BasicExpression::SmithStandardFormRecursive(int defvarindex,
						 int addedvarstartindex,
						 string defvarname,
						 vector<int>& oplabels,
						 vector<Expression>& schemata,
						 vector<Expression>& defcons) {
  // standardize depth-first
  int sz = GetSize();
  if (!IsLeaf()) {
    int defconsize;
    for(int i = 0; i < sz; i++) {
      defconsize = defcons.size();
      GetNode(i)->SmithStandardFormRecursive(defvarindex, addedvarstartindex,
					     defvarname, oplabels, schemata,
					     defcons);
      defvarindex += defcons.size() - defconsize;
    }
  }
  // do this node
  Expression e(0.0);
  bool done = false;
  int dvi = 0;
  double theCoeff;
  for(vector<Expression>::iterator schi = schemata.begin();
      schi != schemata.end() && !done; schi++) {
    if (IsEqualBySchema(*schi)) {
      done = true;
      // save coefficient of root node
      theCoeff = GetCoeff();
      // and re-set it to one
      SetCoeff(1.0);
      // see if *this is already in defcons
      dvi = SmithFindInDefcons(addedvarstartindex, defcons);
      if (dvi == -1) {
	// it isn't, save this into defcons
	e.SetTo(*this);
	defcons.push_back(e);
      } else {
	defvarindex = dvi;
      }
      // replace this node with w_defvarindex
      // -- segvs (?)
      for(int i = 0; i < GetSize(); i++) {
	RecursiveDestroy(GetNodePtr(i));
      }
      DeleteAllNodes();
      SetOpType(VAR);
      SetVarIndex(defvarindex);
      std::stringstream ss;
      ss << defvarname << defvarindex;
      SetVarName(ss.str());
      SetCoeff(theCoeff);
      SetExponent(1.0);
      break;
    }
  }
  if (!done) {
    for(vector<int>::iterator opi = oplabels.begin();
	opi != oplabels.end() && !done; opi++) {
      if (IsEqualByOperator(*opi)) {
	done = true;
	// save coefficient of root node
	theCoeff = GetCoeff();
	// and re-set it to one
	SetCoeff(1.0);
	if (*opi == PRODUCT && sz > 2 && sz!=4) {
	  // replace above if with below only if not using RQuarticConvex
        //if (*opi == PRODUCT && sz > 2) {
	  // treat this case separately -- only case of n-ary product
	  // whose Smith reformulation is binary: w = x1 (x2 ... (x_{n-1}x_n))
	  e = GetNode(sz - 1) * GetNode(sz - 2);
	  // see if e is already in defcons
	  dvi = e->SmithFindInDefcons(addedvarstartindex, defcons);
	  if (dvi == -1) {
	    // it isn't, save
	    defcons.push_back(e);
	  } else {
	    defvarindex = dvi;
	  }
	  for(int i = sz-3; i >= 0; i--) {
	    std::stringstream ss;
	    ss << defvarname << defvarindex;
	    Expression wi(1.0, defvarindex, ss.str());
	    Expression e2 = GetNode(i) * wi;
	    // see if already present in defcons
	    dvi = e->SmithFindInDefcons(addedvarstartindex, defcons);
	    if (dvi == -1) {
	      // it isn't, save
	      defcons.push_back(e2);
	      defvarindex++;
	    } else {
	      defvarindex = dvi;
	    }
	  }
	} else {
	  // see if e is already in defcons
	  dvi = SmithFindInDefcons(addedvarstartindex, defcons);
	  if (dvi == -1) {
	    // it isn't, save
	    e.SetTo(*this);
	    defcons.push_back(e);
	  } else {
	    defvarindex = dvi;
	  }
	}
	// replace this node with w_defvarindex
	/* // -- segvs (why does this and not the other similar ones?)
	for(int i = 0; i < GetSize(); i++) {
	  RecursiveDestroy(GetNodePtr(i));
	}
	*/
	DeleteAllNodes();
	SetOpType(VAR);
	SetVarIndex(defvarindex);
	std::stringstream ss;
	ss << defvarname << defvarindex;
	SetVarName(ss.str());
	SetCoeff(theCoeff);
	SetExponent(1.0);
	break;
      }
    }
  }
}

// defvarname is "w", var names are cat(defvarname, vi)
int BasicExpression::ProdBinCont(int defvarindex, string defvarname,
				 map<int,bool>& integrality,
				 map<int,double>& vlb,
				 map<int,double>& vub,
	  		         map<int,string>& vstr,
				 map<int,pair<double,double> >& addvarbounds,
				 vector<Expression>& defcons) {
  int ret = 0;
  addvarbounds.erase(addvarbounds.begin(), addvarbounds.end());
  ret = ProdBinContRecursive(defvarindex, defvarname, integrality,
			     vlb, vub, vstr, addvarbounds, defcons);
  return ret;
}

// defvarname is "w", var names are cat(defvarname, vi)
int BasicExpression::ProdBinContRecursive(int defvarindex, string defvarname,
					  map<int,bool>& integrality,
					  map<int,double>& vlb,
					  map<int,double>& vub,
  		                          map<int,string>& vstr,
					  map<int,pair<double,double> >&
					  addvarbounds,
					  vector<Expression>& defcons){
  int ret = 0;
  // apply depth-first
  int sz = GetSize();
  if (!IsLeaf()) {
    int defconsize;
    for(int i = 0; i < sz; i++) {
      defconsize = defcons.size();
    ret += GetNode(i)->ProdBinContRecursive(defvarindex, defvarname,
					      integrality, vlb, vub, vstr,
					      addvarbounds, defcons);
    }
    defvarindex += ret;
  }
  // do this node
  if (GetOpType() == PRODUCT && GetSize() == 2 &&
      GetNode(0)->GetOpType() == VAR && GetNode(1)->GetOpType() == VAR) {
    // x*y, check binary
    int v1 = GetNode(0)->GetVarIndex();
    int v2 = GetNode(1)->GetVarIndex();
    if ((integrality[v1] && vlb[v1] == 0 && vub[v1] == 1) ||
	(integrality[v2] && vlb[v2] == 0 && vub[v2] == 1)) {
      // one of the vars is binary, call it x and call the other y
      int xi = (integrality[v1] && vlb[v1] == 0 && vub[v1] == 1) ? v1 : v2;
      int yi = (xi == v1) ? v2 : v1;
      // compute bounds for added variable
      pair<double,double> bnd(vlb[yi], vub[yi]);
      bool isfound = false;
      for(map<int,pair<double,double> >::iterator mi = addvarbounds.begin();
	  mi != addvarbounds.end(); mi++) {
	if (mi->first == defvarindex) {
	  isfound = true;
	  break;
	}
      }
      if (isfound) {
	pair<double,double> oldbnd = addvarbounds[defvarindex];
	if (oldbnd.first > bnd.first) {
	  bnd.first = oldbnd.first;
	}
	if (oldbnd.second < bnd.second) {
	  bnd.second = oldbnd.second;
	}
      }
      addvarbounds[defvarindex] = bnd;
      // compute linearization expressions
      std::stringstream ss;
      ss << defvarname << defvarindex;
      Expression w(1.0, defvarindex, ss.str());
      Expression U(vub[yi]);
      Expression L(vlb[yi]);
      Expression x(1.0, xi, vstr[xi]);
      Expression y(1.0, yi, vstr[yi]);
      Expression xu= U * x;
      Expression xl= L * x;
      Expression c1 = w - xu;
      Expression c2 = xl - w;
      Expression c3 =  w - y + L - xl;
      Expression c4 = y -w - U + xu;
      Simplify(&c1);
      Simplify(&c2);
      Simplify(&c3);
      Simplify(&c4);
      defcons.push_back(c1);
      defcons.push_back(c2);
      defcons.push_back(c3);
      defcons.push_back(c4);
      ret++;
      // now replace this product with w
      // -- segvs (?)
      for(int i = 0; i < GetSize(); i++) {
	RecursiveDestroy(GetNodePtr(i));
      }
      DeleteAllNodes();
      SetOpType(VAR);
      SetVarIndex(defvarindex);
      std::stringstream ss2;
      ss2 << defvarname << defvarindex;
      SetVarName(ss2.str());
      SetCoeff(1.0);
      SetExponent(1.0);
    }
  }
  return ret;
}

// vn is "w", var names are cat(vn, vi)
void BasicExpression::ResetVarNames(string vn, int lid, int uid) {
  // set all variable names in the expression to vn
  if (!IsLeaf()) {
    for(int i = 0; i < GetSize(); i++) {
      GetNode(i)->ResetVarNames(vn, lid, uid);
    }
  } else {
    if (GetOpType() == VAR) {
      int vi = GetVarIndex();
      if (vi >= lid && vi <= uid) {
	std::stringstream ss;
	ss << vn << vi;
	SetVarName(ss.str());
      }
    }
  }
}


bool BasicExpression::DistributeProductsOverSums(void) {
  // recursive part
  bool ret = false;
  if (!IsLeaf()) {
    for(int i = 0; i < GetSize(); i++) {
      bool haschanged = GetNode(i)->DistributeProductsOverSums();
      if (haschanged) {
	ret = true;
      }
    }
  }
  // deal with this node
  Expression e(0.0);
  if (GetOpType() == PRODUCT) {
    for(int i = 0; i < GetSize(); i++) {
      if (GetNode(i)->GetOpType() == SUM) {
	// found occurrence of *(+), distribute
	ret = true;
	Expression f = (*this) / GetNode(i);
	Simplify(&f);
	for(int j = 0; j < GetNode(i)->GetSize(); j++) {
	  e = e + f * GetNode(i)->GetNode(j);
	}
	// now replace this with e
	ReplaceWithExpression(e);
      }
    }
  }
  return ret;
}

void BasicExpression::GetVarIndices(vector<int>& vidx) {
  if (!IsLeaf()) {
    for(int i = 0; i < GetSize(); i++) {
      GetNode(i)->GetVarIndices(vidx);
    }
  } else if (IsVariable()) {
    int vi = GetVarIndex();
    if (find(vidx.begin(), vidx.end(), vi) == vidx.end()) {
      vidx.push_back(vi);
    }
  }
}

void BasicExpression::GetVarIndicesInSchema(vector<int>& vidx,
					    Expression schema) {
  // recurse
  if (!IsLeaf()) {
    for(int i = 0; i < GetSize(); i++) {
      GetNode(i)->GetVarIndicesInSchema(vidx, schema);
    }
  }
  // deal with this node
  if (IsEqualBySchema(schema)) {
    GetVarIndices(vidx);
  }
}

// find the variable name corresponding to variable index vi
string BasicExpression::FindVariableName(int vi) {
  string vn;
  if (IsVariable()) {
    if (GetVarIndex() == vi) {
      return GetVarName();
    } else {
      return "";
    }
  } else {
    int i;
    for(i = 0; i < GetSize(); i++) {
      vn = GetNode(i)->FindVariableName(vi);
      if (vn.length() > 0) {
	return vn;
      }
    }
  }
  return "";
}

// is this expression linear?
bool BasicExpression::IsLinear(void) const {
  if (IsVariable()) {
    if (GetExponent() != 0 && GetExponent() != 1) {
      return false;
    } else {
      return true;
    }
  }
  if (IsConstant()) {
    return true;
  }
  if(GetOpType() == SUM || GetOpType() == DIFFERENCE) {
    int i;
    for(i = 0; i < GetSize(); i++) {
      if (!GetNode(i)->IsLinear()) {
	return false;
      }
    }
    return true;
  } else {
    return false;
  }
}

// is this expression a quadratic product?
bool BasicExpression::IsQuadratic(int& prodtype) const {
  bool ret = false;
  if ((GetOpType() == PRODUCT &&
       GetNode(0)->GetOpType() == VAR && GetNode(1)->GetOpType() == VAR) ||
      (GetOpType() == POWER &&
       GetNode(0)->GetOpType() == VAR && GetNode(1)->GetValue() == 2) ||
      (GetOpType() == VAR && GetExponent() == 2)) {
    prodtype = GetOpType();
    ret = true;
  }
  return ret;
}
bool BasicExpression::IsQuadratic(void) const {
  int ret = 0;
  return IsQuadratic(ret);
}

// return info about the linear part (assumes Simplify() has already been
// called on this) - return false if expression has no linear part
// by "linear part" we mean lin(x) in expr(x,y) = lin(x) + nonlin(y)
bool BasicExpression::GetLinearInfo(vector<double>& lincoeff,
				    vector<int>& linvi,
				    vector<string>& linvn,
				    double& c) {
  c = 0;
  bool ret = false;
  Expression nl(GetPureNonlinearPart());
  if (IsLinear()) {
    nl->Zero();
  }
  if (lincoeff.size() > 0) {
    lincoeff.erase(lincoeff.begin(), lincoeff.end());
    linvi.erase(linvi.begin(), linvi.end());
    linvn.erase(linvn.begin(), linvn.end());
  }
  if (IsLeaf()) {
    if (IsConstant()) {
      // just a constant
      c = GetValue();
      ret = true;
    } else if (IsVariable() && GetExponent() == 1) {
      // just a variable
      linvi.push_back(GetVarIndex());
      lincoeff.push_back(GetCoeff());
      linvn.push_back(GetVarName());
      ret = true;
    }
  } else {
    if (GetOpType() == SUM) {
      int i, vi;
      c = 0;
      for(i = 0; i < GetSize(); i++) {
	if (GetNode(i)->IsConstant()) {
	  if (i > 0) {
	    cerr << "BasicExpression::GetLinearInfo: WARNING: "
		 << "run Simplify() first\n";
	  }
	  c += GetNode(i)->GetValue();
	} else if (GetNode(i)->IsVariable() &&
		   GetNode(i)->GetExponent() == 1) {
	  vi = GetNode(i)->GetVarIndex();
	  if (!nl->DependsOnVariable(vi)) {
	    // vi depends only linearly on expression
	    linvi.push_back(vi);
	    lincoeff.push_back(GetNode(i)->GetCoeff());
	    linvn.push_back(GetNode(i)->GetVarName());
	    ret = true;
	  }
	}
      }
    }
  }
  return ret;
}

// return info about the purely linear part
// (assumes Simplify() has already been
// called on this) - return false if expression has no linear part
// by "pure linear part" we mean e.g. x+y in x+y+y^2
bool BasicExpression::GetPureLinearInfo(vector<double>& lincoeff,
					vector<int>& linvi,
					vector<string>& linvn,
					double& c) {
  c = 0;
  bool ret = false;
  int i;
  Expression nl(GetPureNonlinearPart());
  if (IsLinear()) {
    nl->Zero();
  }
  if (lincoeff.size() > 0) {
    lincoeff.erase(lincoeff.begin(), lincoeff.end());
    linvi.erase(linvi.begin(), linvi.end());
    linvn.erase(linvn.begin(), linvn.end());
  }
  if (IsLeaf()) {
    if (IsConstant()) {
      // just a constant
      c = GetValue();
      ret = true;
    } else if (IsVariable() && GetExponent() == 1) {
      // just a variable
      linvi.push_back(GetVarIndex());
      lincoeff.push_back(GetCoeff());
      linvn.push_back(GetVarName());
      ret = true;
    }
  } else {
    if (GetOpType() == SUM) {
      int vi;
      c = 0;
      for(i = 0; i < GetSize(); i++) {
	if (GetNode(i)->IsConstant()) {
	  if (i > 0) {
	    cerr << "BasicExpression::GetLinearInfo: WARNING: "
		 << "run Simplify() first\n";
	  }
	  c += GetNode(i)->GetValue();
	} else if (GetNode(i)->IsVariable() &&
		   GetNode(i)->GetExponent() == 1) {
	  vi = GetNode(i)->GetVarIndex();
	  linvi.push_back(vi);
	  lincoeff.push_back(GetNode(i)->GetCoeff());
	  linvn.push_back(GetNode(i)->GetVarName());
	  ret = true;
	}
      }

    }
  }
  return ret;
}

// get the linear part - x in x+y+y^2
Expression BasicExpression::GetLinearPart(void) {
  vector<double> lincoeff;
  vector<int> linvi;
  vector<string> linvn;
  double c;
  GetLinearInfo(lincoeff, linvi, linvn, c);
  int i;
  Expression ret;
  if (lincoeff.size() > 0) {
    ret->SetOpType(VAR);
    ret->SetVarIndex(linvi[0]);
    ret->SetCoeff(lincoeff[0]);
    ret->SetVarName(linvn[0]);
    ret->SetExponent(1);
    if (lincoeff.size() > 1) {
      Expression addend(1.0, -1, NOTVARNAME);
      for(i = 1; i < (int) lincoeff.size(); i++) {
	addend->SetVarIndex(linvi[i]);
	addend->SetCoeff(lincoeff[i]);
	addend->SetVarName(linvn[i]);
	ret = ret + addend;
      }
    }
  }
  return ret;
}

// get the pure linar part - x+y in x+y+y^2
Expression BasicExpression::GetPureLinearPart(void) {
  vector<double> lincoeff;
  vector<int> linvi;
  vector<string> linvn;
  double c;
  GetPureLinearInfo(lincoeff, linvi, linvn, c);
  int i;
  Expression ret(0.0);
  if (lincoeff.size() > 0) {
    ret->SetOpType(VAR);
    ret->SetVarIndex(linvi[0]);
    ret->SetCoeff(lincoeff[0]);
    ret->SetVarName(linvn[0]);
    ret->SetExponent(1);
    if (lincoeff.size() > 1) {
      for(i = 1; i < (int) lincoeff.size(); i++) {
	Expression addend(lincoeff[i], linvi[i], linvn[i]);
	ret = SumLink(ret, addend);
      }
    }
  }
  return ret;
}

// get the nonlinear part - nonlin(y) in expr(x,y) = lin(x) + nonlin(y)
Expression BasicExpression::GetNonlinearPart(void) {
  Expression ret(GetPureNonlinearPart());
  vector<double> linval;
  vector<int> linidx;
  vector<string> linvn;
  double c = 0;
  GetPureLinearInfo(linval, linidx, linvn, c);
  int i;
  Expression addend(1.0, -1, NOTVARNAME);
  // we cycle backwards to keep the order of addends in ret
  for(i = linidx.size() - 1; i >= 0 ; i--) {
    if (ret->DependsOnVariable(linidx[i])) {
      // pure nonlinear part depends on varindex i, add it to ret
      addend->SetCoeff(linval[i]);
      addend->SetVarIndex(linidx[i]);
      addend->SetVarName(linvn[i]);
      ret = addend + ret;
    }
  }
  return ret;
}


// get the purely nonlinear part - e.g. only y^2 in x+y+y^2
Expression BasicExpression::GetPureNonlinearPart(void) {
  Expression ret(0.0);
  if (!IsLeaf()) {
    if (GetOpType() == SUM) {
      int i;
      for(i = 0; i < GetSize(); i++) {
	if (!GetNode(i)->IsLinear()) {
	  ret = SumLink(ret, GetNode(i));
	}
      }
    } else if (GetOpType() == DIFFERENCE) {
      int i;
      for(i = 0; i < GetSize(); i++) {
	if (!GetNode(i)->IsLinear()) {
	  ret = ret - GetNode(i);
	}
      }
    } else if (GetOpType() == PLUS) {
      ret = GetNode(0);
    } else if (GetOpType() == MINUS) {
      ret = GetNode(0);
      ret->SetCoeff(- ret->GetCoeff());
    } else {
      ret = *this;
    }
  } else {
    // leaf but can be a power
    if (GetExponent() != 0 && GetExponent() != 1) {
      ret = *this;
    }
  }
  return ret;
}

// get value of additive constant
double BasicExpression::GetConstantPart(void) {
  double ret = 0;
  if (IsConstant()) {
    ret = GetValue();
  } else if (!IsLeaf()) {
    int op = GetOpType();
    if (op == SUM || op == DIFFERENCE) {
      int i = 0;
      int sz = GetSize();
      while(i < sz) {
	if (GetNode(i)->IsConstant()) {
	  if (op == SUM || (op == DIFFERENCE && i == 0)) {
	    ret += GetNode(i)->GetValue();
	  } else {
	    ret -= GetNode(i)->GetValue();
	  }
	}
	i++;
      }
    }
  }
  return ret;
}

// doesn't deal with additive constants in PLUS/MINUS operands
double BasicExpression::RemoveAdditiveConstant(void) {
  double ret = 0;
  if (IsConstant()) {
    ret = GetValue();
    SetValue(0);
  } else if (!IsLeaf()) {
    int op = GetOpType();
    if (op == SUM || op == DIFFERENCE) {
      int i = 0;
      int sz = GetSize();
      while(i < sz) {
	if (GetNode(i)->IsConstant()) {
	  if (op == SUM || (op == DIFFERENCE && i == 0)) {
	    ret += GetNode(i)->GetValue();
	  } else {
	    ret -= GetNode(i)->GetValue();
	  }
	  DeleteNode(i);
	  sz--;
	} else {
	  i++;
	}
      }
    }
  }
  return ret;
}

// perform interval arithmetics on the expression;
// Vlb[i] is the lower bound of varindex i, Vub is similar,
// returned elb and eub hold values to expression interval
void BasicExpression::Interval(map<int,double> Vlb, map<int,double> Vub,
			       double& elb, double& eub) {
  int i;
  if (IsLeaf()) {
    // leaf
    if (IsConstant()) {
      elb = GetValue();
      eub = elb;
    } else if (IsVariable()) {
      double expon = GetExponent();
      double coeff = GetCoeff();
      int vi = GetVarIndex();
      if (expon == 1) {
	// exponent is 1, simple variable
	elb = coeff * Vlb[vi];
	eub = coeff * Vub[vi];
      } else {
	// var ^ constant
	constpowermkrange(Vlb[vi], Vub[vi], expon, &elb, &eub);
	elb = coeff * elb;
	eub = coeff * eub;
      }
    }
  } else {
    // operator recurse
    vector<double> tmplb, tmpub;
    double tlb, tub;
    for(i = 0; i < GetSize(); i++) {
      tlb = -Ev3Infinity();
      tub = Ev3Infinity();
      GetNode(i)->Interval(Vlb, Vub, tlb, tub);
      tmplb.push_back(tlb);
      tmpub.push_back(tub);
    }
    //// changed 080722 to take into account a possible coefficient
    double coeff = GetCoeff();
    double cf = 1.;
    if (coeff != 1) cf = cf*coeff;
    if(cf != 1.) {
      bilinearprodmkrange(tmplb[0], tmpub[0], cf, cf, &tmplb[0], &tmpub[0]);
    }
    ////
    // do the gig
    int op = GetOpType();
    int sz = GetSize();
    double t1 = 0;
    double t2 = 0;
    double t3 = 0;
    double t4 = 0;
    double t5 = 0;
    switch(op) {
    case SUM:
      for(i = 0; i < sz; i++) {
	t1 += tmplb[i];
	t2 += tmpub[i];
      }
      break;
    case DIFFERENCE:
      for(i = 1; i < sz; i++) {
	t1 += tmplb[i];
	t2 += tmpub[i];
      }
      t1 = tmplb[0] - t2;
      t2 = tmpub[0] - t1;
      break;
    case PRODUCT:
      assert(sz > 1);
      bilinearprodmkrange(tmplb[0], tmpub[0], tmplb[1], tmpub[1], &t1, &t2);
      for(i = 2; i < sz; i++) {
	t3 = t1;
	t4 = t2;
	bilinearprodmkrange(t3, t4, tmplb[i], tmpub[i], &t1, &t2);
      }
      break;
    case FRACTION:
      fractionmkrange(tmplb[0], tmpub[0], tmplb[1], tmpub[1], &t1, &t2);
      break;
    case POWER:
      powermkrange(tmplb[0], tmpub[0], tmplb[1], tmpub[1], &t1, &t2);
      break;
    case MINUS:
      t1 = -tmpub[0];
      t2 = -tmplb[0];
      break;
    case LOG:
      if (tmplb[0] <= 0 && tmpub[0] <= 0) {
	t1 = 0;
	t2 = 0;
      } else if (tmplb[0] <= 0) {
	t1 = -Ev3Infinity();
	t2 = std::log(tmpub[0]);
      } else {
	t1 = std::log(tmplb[0]);
	t2 = std::log(tmpub[0]);
      }
      break;
    case EXP:
      t1 = std::exp(tmplb[0]);
      t2 = std::exp(tmpub[0]);
      break;
    case SIN:
      t3 = tmpub[0] - tmplb[0];
      t4 = tmplb[0] / PEV3PI + 0.5;
      t5 = tmpub[0] / PEV3PI + 0.5;
      if (t3 < PEV3PI && std::floor(t4) == std::floor(t5)) {
	t1 = std::sin(tmplb[0]);
	t2 = std::sin(tmpub[0]);
	if (t1 > t2) {
	  double ttmp = t1;
	  t1 = t2;
	  t2 = ttmp;
	}
      } else {
	// should provide an "else if" for case where we have
	// -1 < v < UB or LB < v < 1 but I can't be bothered
	t1 = -1;
	t2 = 1;
      }
      break;
    case COS:
      t3 = tmpub[0] - tmplb[0];
      t4 = tmplb[0] / PEV3PI + 0.5;
      t5 = tmpub[0] / PEV3PI + 0.5;
      if (t3 < PEV3PI && std::floor(t4) == std::floor(t5)) {
	t1 = std::cos(tmplb[0]);
	t2 = std::cos(tmpub[0]);
	if (t1 > t2) {
	  double ttmp = t1;
	  t1 = t2;
	  t2 = ttmp;
	}
      } else {
	// should provide an "else if" for case where we have
	// -1 < v < UB or LB < v < 1 but I can't be bothered
	t1 = -1;
	t2 = 1;
      }
      break;
    case TAN:
      t3 = tmpub[0] - tmplb[0];
      t4 = tmplb[0] / PEV3PI + 0.5;
      t5 = tmpub[0] / PEV3PI + 0.5;
      if (t3 < PEV3PI && std::floor(t4) == std::floor(t5)) {
	// lb and ub within the same tangent branch
	t1 = std::tan(tmplb[0]);
	t2 = std::tan(tmpub[0]);
	if (t1 > t2) {
	  double ttmp = t1;
	  t1 = t2;
	  t2 = ttmp;
	}
      } else {
	// should provide an "else if" for case where we have
	// -1 < v < UB or LB < v < 1 but I can't be bothered
	t1 = -Ev3Infinity();
	t2 = Ev3Infinity();
      }
      break;
    case COT:
      cerr << "expression.cxx::Interval(): cot() not implemented\n";
      break;
    case SINH:
      cerr << "expression.cxx::Interval(): sinh() not implemented\n";
      break;
    case COSH:
      cerr << "expression.cxx::Interval(): cosh() not implemented\n";
      break;
    case TANH:
      cerr << "expression.cxx::Interval(): tanh() not implemented\n";
      break;
    case COTH:
      cerr << "expression.cxx::Interval(): coth() not implemented\n";
      break;
    case SQRT:
      t1 = std::sqrt(tmplb[0]);
      t2 = std::sqrt(tmpub[0]);
      break;
    default:
      break;
    }
    elb = t1;
    eub = t2;
  }
  if (elb > eub) {
    // switch if needed (it might happen in certain cases that we end
    // up with an inverted range)
    double tmp = elb;
    elb = eub;
    eub = tmp;
  }
  SetLB(elb);
  SetUB(eub);
}

void BasicExpression::FBBTUpDown(map<int,double>& Vlb, map<int,double>& Vub, 
				 double elb, double eub) {
  double ielb = -LARGE;
  double ieub = LARGE;
  // UP
  Interval(Vlb, Vub, ielb, ieub);
  if (ielb >= elb && ieub <= eub) {
    // FBBT has no effecton this expression
    return;
  }
  // [elb,eub] = UP \cap [elb,eub]
  if (ielb > elb) {
    elb = ielb;
  }
  if (ieub < eub) {
    eub = ieub;
  }
  // DOWN
  FBBTDown(Vlb, Vub, elb, eub);
}

void BasicExpression::FBBTDown(map<int,double>& Vlb, map<int,double>& Vub,
			       double elb, double eub) {
  int k;
  double l0, u0, l1, u1, l2, u2;
  double c;
  if (!IsLeaf()) {
    // can only do DOWN on non-leaf nodes
    k = GetSize();
    // case DIFFERENCE 
    if (GetOpType() == DIFFERENCE && k == 2) {
      // inverse operator is SUM (of two operands)
      l0 = GetNode(0)->GetLB();
      u0 = GetNode(0)->GetUB();
      l1 = GetNode(1)->GetLB();
      u1 = GetNode(1)->GetUB();
      // set bounds for left node 
      GetNode(0)->SetLB(max(l0, elb + l1));
      GetNode(0)->SetUB(min(u0, eub + u1));
      // set bounds for right node
      GetNode(1)->SetLB(max(l1, elb + l0));
      GetNode(1)->SetUB(min(u1, eub + u0));
    } else if (GetOpType() == SUM && k >= 2) {
      // inverse operator is DIFFERENCE (top - sum(subnodes but one))
      // compute sum of all intervals in [l0,u0] 
      l0 = 0;
      u0 = 0;
      for(int i = 0; i < k; i++) {
	l0 += GetNode(i)->GetLB();
	u0 += GetNode(i)->GetUB();
      }
      // assign new intervals to each node
      for(int i = 0; i < k; i++) {
	l1 = GetNode(i)->GetLB();
	u1 = GetNode(i)->GetUB();
	// now perform top node interval minus the sum of all interv but i-th
	// ([l0-l1,u0-u1] is the sum of all intervals but the i-th)
	l2 = max(l1, elb - (u0 - u1));
	u2 = min(u1, eub - (l0 - l1));
	GetNode(i)->SetLB(l2);
	GetNode(i)->SetUB(u2);
      }
    }
    // recurse
    for(int i = 0; i < k; i++) {
      GetNode(i)->FBBTDown(Vlb, Vub, GetNode(i)->GetLB(), GetNode(i)->GetUB());
    }
  } else {
    // copy node ranges into variable ranges
    if (IsVariable()) {
      k = GetVarIndex();
      c = GetExponent();
      l1 = GetCoeff();
      if (l1 > 0) {
	l0 = GetLB() / l1;
        u0 = GetUB() / l1;
      } else if (l1 < 0) {
	l0 = GetUB() / l1;
	u0 = GetLB() / l1;
      } else {
	l0 = 0;
	u0 = 0;
      }
      if (l1 != 0) {
	int ci = (int) rint(c);
	if (fabs(c - ci) < 1/LARGE && c >= 0) {
	  // nonnegative integer exponent
	  if (ci == 0) {
	    // zero
	    l0 = 1;
	    u0 = 1;
	  } else if (ci == 1) {
	    // stay with current values for l0, u0
	  } else if (ci % 2 == 0) {
	    // even not in {0,1}
	    if (l0 <= 0) {
	      l0 = 0;
	    } else {
	      l0 = pow(l0, 1/c);
	    }
	    if (u0 <= 0) {
	      u0 = 0;
	    } else {
	      u0 = pow(u0, 1/c);
	    }
	  } else {
	    // odd not in {0,1}
	    l0 = pow(l0, 1/c);
	    u0 = pow(u0, 1/c);
	  }
	} else {
	  cerr << "FBBTDown(): cannot deal with negative exponents" << endl;
	  exit(65);
	}
      }
      if (l0 > Vlb[k]) {
	Vlb[k] = l0; 
      }
      if (u0 < Vub[k]) {
	Vub[k] = u0;
      }
    }
  }
}

int BasicExpression::TreeData(int r, 
			      std::vector<std::pair<int,int> >& edgeList, 
			      std::map<int,int>& nodeType, 
			      std::map<int,int>& index, 
			      std::map<int,double>& coeff,
			      std::map<int,pair<double,double> >& bnds) {
  edgeList.erase(edgeList.begin(), edgeList.end());
  nodeType.erase(nodeType.begin(), nodeType.end());
  index.erase(index.begin(), index.end());
  coeff.erase(coeff.begin(), coeff.end());
  bnds.erase(bnds.begin(), bnds.end());
  int ret = 1;
  TreeRecursive(1, ret, edgeList, nodeType, index, coeff, bnds);
  return ret;
}

void BasicExpression::TreeRecursive(int r, int& maxr,
				    std::vector<std::pair<int,int> >& edgeList, 
				    std::map<int,int>& nodeType, 
				    std::map<int,int>& index, 
				    std::map<int,double>& coeff,
				    std::map<int,std::pair<double,double> >& 
				    bnds) {
  pair<int,int> p(r, maxr);
  if (r != maxr) {
    edgeList.push_back(p);
  }

  if (IsLeaf()) {
    if (IsVariable()) {
      // variable leaf
      double e = GetExponent();
      double c = GetCoeff();
      int vi = GetVarIndex();
      if (fabs(e-1) < 1/LARGE && fabs(c-1) < 1/LARGE) {
	// exponent = coeff = 1, make node 1*var
	// add product node
	nodeType[maxr] = PRODUCT;
	p.first = maxr;
	bnds[maxr].first = GetLB();
	bnds[maxr].second = GetUB();
	// add coeff node
	maxr++;
	nodeType[maxr] = CONST;
	coeff[maxr] = 1;
	bnds[maxr].first = 1;
	bnds[maxr].second = 1;
	p.second = maxr;
	edgeList.push_back(p);
	// add var node
	maxr++;
	nodeType[maxr] = VAR;
	index[maxr] = vi;
	bnds[maxr].first = GetLB();
	bnds[maxr].second = GetUB();
	p.second = maxr;
	edgeList.push_back(p);
      } else if (fabs(e-1) < 1/LARGE) {
	// exponent = 1 but coeff != 1
	// add product node
	nodeType[maxr] = PRODUCT;
	p.first = maxr;
	if (c >= 0) {
	  bnds[maxr].first = GetLB() * c;
	  bnds[maxr].second = GetUB() * c;
	} else {	
	  bnds[maxr].first = GetUB() * c;
	  bnds[maxr].second = GetLB() * c;
	}  
	// add coeff node
	maxr++;
	nodeType[maxr] = CONST;
	coeff[maxr] = c;
	bnds[maxr].first = c;
	bnds[maxr].second = c;
	p.second = maxr;
	edgeList.push_back(p);
	// add var node
	maxr++;
	nodeType[maxr] = VAR;
	index[maxr] = vi;
	bnds[maxr].first = GetLB();
	bnds[maxr].second = GetUB();
	p.second = maxr;
	edgeList.push_back(p);
      } else if (fabs(c - 1) < 1/LARGE) {
	// coeff = 1 but exponent != 1
	// add power node
	nodeType[maxr] = POWER;
	p.first = maxr;
	cout << "TreeRecursive(): fake UP not implemented for exponents\n";
	exit(69);
	// add var node
	maxr++;
	nodeType[maxr] = VAR;
	index[maxr] = vi;
	bnds[maxr].first = GetLB();
	bnds[maxr].second = GetUB();
	p.second = maxr;
	edgeList.push_back(p);
	// add exponent node
	maxr++;
	nodeType[maxr] = CONST;
	coeff[maxr] = e;
	bnds[maxr].first = e;
	bnds[maxr].second = e;
	p.second = maxr;
	edgeList.push_back(p);
      } else {
	// coeff and exponent != 1
	// add product node
	nodeType[maxr] = PRODUCT;
	p.first = maxr;
	cout << "TreeRecursive(): fake UP not implemented for exponents\n";
	exit(69);
	// add coeff node
	maxr++;
	nodeType[maxr] = CONST;
	coeff[maxr] = c;
	bnds[maxr].first = c;
	bnds[maxr].second = c;
	p.second = maxr;
	edgeList.push_back(p);
	// add power node
	maxr++;
	nodeType[maxr] = POWER;
	cout << "TreeRecursive(): fake UP not implemented for exponents\n";
	exit(69);
	p.first = maxr;
	// add var node
	maxr++;
	nodeType[maxr] = VAR;
	index[maxr] = vi;
	bnds[maxr].first = GetLB();
	bnds[maxr].second = GetUB();
	p.second = maxr;
	edgeList.push_back(p);
	// add exponent node
	maxr++;
	nodeType[maxr] = CONST;
	coeff[maxr] = e;
	bnds[maxr].first = e;
	bnds[maxr].second = e;
	p.second = maxr;
	edgeList.push_back(p);
      }
    } else {
      // constant leaf
      nodeType[maxr] = GetOpType();
      coeff[maxr] = GetValue();
    }

  } else {
    // recursive part
    nodeType[maxr] = GetOpType();
    for(int i = 0; i < GetSize(); i++) {
      maxr++;
      GetNode(i)->TreeRecursive(r, maxr, edgeList, 
				nodeType, index, coeff, bnds);
    }    
  }
}

bool BasicExpression::IsSmithStandard(void) {
  return IsSmithStandard(0, 0);
}

bool BasicExpression::IsSmithStandard(double l, double u) {
  bool ret = false;
  if (!IsLeaf()) {
    int op = GetOpType();
    if (GetSize() == 2 && (op == SUM || op == DIFFERENCE)) {
      if (!GetNode(0)->IsConstant() && !GetNode(1)->IsConstant()) {
	int v1 = -1;
	int v2 = -1;
	if (GetNode(0)->IsVariable()) {
	  v1 = GetNode(0)->GetVarIndex();
	}
	if (GetNode(1)->IsVariable()) {
	  v2 = GetNode(1)->GetVarIndex();
	}
	if (!(v1 > -1 && v2 > -1)) {
	  int posv = -1;
	  int posop = -1;
	  if (v1 > -1) {
	    posv = 0;
	    posop = 1;
	  } else if (v2 > -1) {
	    posv = 1;
	    posop = 0;
	  }
	  if (posv > -1) {
	    if (GetNode(posv)->GetExponent() == 1) {
	      // we have finally determined that expression is of the form
	      // c*w_i +/- something or something +/- c*w_i, that
	      // something is an expression which is not a leaf, and that
	      // it is in GetNode(pos)
	      Expression e = GetNode(posop);
#ifdef NOSMITH1VARTRIPLES
	      int vi = GetNode(posv)->GetVarIndex();
	      if (e->DependsOnVariable(vi)) {
		// the nonlinear part depends on the "added variable",
		// this is not a good "triple", must reformulate
		ret = false;
	      } else
#endif
		{
		if (!(e->IsLinear())) {
		  int sz = e->GetSize();
		  if (sz == 1 && e->GetNode(0)->IsLeaf()) {
		    // unary function with leaf argument: OK
		    ret = true;
		  } else if (sz == 2 && e->GetNode(0)->IsLeaf() &&
			     e->GetNode(1)->IsLeaf()){
		    // binary operation with two leaves arguments
		    if (e->GetOpType() == PRODUCT) {
		      // binary product with two leaves arguments
		      // and whatever constraint bounds: OK
		      ret = true;
		    } else if (e->GetOpType() == FRACTION &&
			       l == u && u == 0) {
		      // binary fraction with two leaves arguments
		      // and constraint bounds l == u == 0: OK
		      ret = true;
		    } else if (e->GetOpType() == POWER && e->GetSize() == 2 &&
			       e->GetNode(1)->IsConstant()) {
		      if (e->GetNode(1)->GetValue() >= 0) {
			// binary power with constant positive exponent and
			// any constraint bounds: OK
			ret = true;
		      } else {
			if (l == u) {
			  // constant negative exponent and
			  // constraint is equation: OK
			  ret = true;
			} else {
			  // constant negative exponent and
			  // inequality constraint: not ok
			  ret = false;
			}
		      }
		    } else if (e->GetOpType() == LOG ||
			       e->GetOpType() == EXP ||
			       e->GetOpType() == SQRT) {
		      // univariate convex/concave with any
		      // constraint bounds: OK
		      // (can do some more complex testing here to
		      // include univariate convex/concave parts of
		      // trigonometric/hyperbolic functions)
		      ret = true;
		    } else {
		      // anything else is false, needs standardization proper
		      ret = false;
		    }
		  }
		}
	      }
	    }
	  }
	} else {
	  // both addends are variables, check that one has exponent 1
	  // and the other != 1 and != 0
	  double expon1 = GetNode(0)->GetExponent();
	  double expon2 = GetNode(1)->GetExponent();
	  if (expon1 != 0 && expon2 != 0) {
	    if (((expon1 == 1 && expon2 != 1) ||
		 (expon1 != 1 && expon2 == 1))
#ifdef NOSMITH1VARTRIPLES
		&& (v1 != v2)
#endif
		) {
	      ret = true;
	    }
	  }
	}
      }
    }
  }
  return ret;
}

bool BasicExpression::IsOptStandard(void) {
  bool ret = false;
  // not implemented yet
  return ret;
}

bool BasicExpression::IsEvidentlyConvex(void) {
  if (IsLinear()) {
    return true;
  }
  bool ret = true;
  int sz = GetSize();
  int op = GetOpType();
  int i;
  double c = 0;
  double expon = 0;
  switch(op) {
  case SUM:
    for(i = 0; i < sz; i++) {
      if (!GetNode(i)->IsLinear()) {
	if (!GetNode(i)->IsEvidentlyConvex()) {
	  ret = false;
	  break;
	}
      }
    }
    // invert if this' coefficient is negative
    c = GetCoeff();
    if (ret) {
      if (c < 0) {
	ret = false;
      }
    }
    else
    {
      ret = true;
      for(i = 0; i < sz; i++) {
        if (!GetNode(i)->IsLinear()) {
   	  if (!GetNode(i)->IsEvidentlyConcave()) {
 	    ret = false;
	    break;
	  }
        }
      }
      // invert if this' coefficient is negative
      c = GetCoeff();
      if (ret) {
        if (c > 0) {
  	  ret = false;
        }
      }
    }

    break;
  case DIFFERENCE:
    if (!GetNode(0)->IsLinear()) {
      if (!GetNode(0)->IsEvidentlyConvex()) {
	ret = false;
	break;
      }
    }
    if (ret) {
      for(i = 1; i < sz; i++) {
	if (!GetNode(i)->IsLinear()) {
	  if (!GetNode(i)->IsEvidentlyConcave()) {
	    ret = false;
	    break;
	  }
	}
      }
    }
    c = GetCoeff();
    if (ret) {
      if (c < 0) {
	ret = false;
      }
    }
    else
    {
      ret = true;
      if (!GetNode(0)->IsLinear()) {
        if (!GetNode(0)->IsEvidentlyConcave()) {
  	  ret = false;
	  break;
        }
      }
      if (ret) {
        for(i = 1; i < sz; i++) {
	  if (!GetNode(i)->IsLinear()) {
	    if (!GetNode(i)->IsEvidentlyConvex()) {
	      ret = false;
	      break;
	    }
	  }
        }
      }
      c = GetCoeff();
      if (ret) {
        if (c > 0) {
	  ret = false;
        }
      }
    }
    break;
  case PRODUCT: case FRACTION: case SIN: case COS: case TAN: case COT:
  case SINH: case COSH: case TANH: case COTH:
    // in general, neither convex nor concave
    ret = false;
    break;
  case LOG: case SQRT: //CLAUDIA
    // depending on coefficient, Concave (c>0) or convex
    c = GetCoeff();
    if (c < 0) {
      ret = true;
    } else {
      ret = false;
    }
    break;
  case EXP: //CLAUDIA
    // depending on coefficient, Convex (c>0) or concave
    c = GetCoeff();
    if (c < 0) {
      ret = true;
    } else {
      ret = false;
    }
    break;
  case POWER:
    if (GetSize() == 2 && GetNode(1)->IsConstant()) {
      c = GetCoeff();
      expon = GetNode(1)->IsConstant();
      if (!is_odd(expon) && expon > 1) {
	// positive non-odd rational exponent > 1
	if (c > 0) {
	  ret = true;
	} else {
	  ret = false;
	}
      } else if (expon > 0 && expon < 1) {
	// positive 0 < exponent < 1
	if (c < 0) {
	  ret = true;
	} else {
	  ret = false;
	}
      } else {
	// all other cases are in general neither concave nor convex
	// (we don't know the range of the variable here)
	ret = false;
      }
    } else {
      ret = false;
    }
    break;
  case VAR:
    if (GetExponent() != 1) {
      // this is like POWER above
      c = GetCoeff();
      expon = GetExponent();
      if (!is_odd(expon) && expon > 1) {
	// positive non-odd rational exponent > 1
	if (c > 0) {
	  ret = true;
	} else {
	  ret = false;
	}
      } else if (expon > 0 && expon <1) {
	// positive 0 < exponent < 1
	if (c < 0) {
	  ret = true;
	} else {
	  ret = false;
	}
      } else {
	// all other cases are in general neither concave nor convex
	// (we don't know the range of the variable here)
	ret = false;
      }
    } else {
      // this should have been caught by the IsLinear() at the beginning
      // anyway
      ret = true;
    }
  }
  return ret;
}

bool BasicExpression::IsEvidentlyConcave(void) {
  if (IsLinear()) {
    return true;
  }
  bool ret = false;
  double c = GetCoeff();
  SetCoeff(-c);
  ret = IsEvidentlyConvex();
  SetCoeff(c);
  return ret;
}

bool BasicExpression::IsEvidentlyConvex(double l, double u) {
  bool ret = true;
  if (IsLinear()) {
    return true;
  }
  if (l > -Ev3Infinity() && u < Ev3Infinity()) {
    // a nonlinear constraint with both bounds finite always
    // defines a nonconvex set
    ret = false;
  } else if (l <= -Ev3Infinity() && u >= Ev3Infinity()) {
    // this is an inactive constraint: always convex
    ret = true;
  } else if (l <= -Ev3Infinity()) {
    // situation is -inf < expr < u. If expr is convex, then
    // we have a convex constraint
    if (IsEvidentlyConvex()) {
      ret = true;
    } else {
      ret = false;
    }
  } else if (u >= Ev3Infinity()) {
    // situation is l < expr < inf. If expr is concave, then
    // we have a convex constraint
    if (IsEvidentlyConcave()) {
      ret = true;
    } else {
      ret = false;
    }
  }
  return ret;
}

// this adds the expression tree to an existing expressions DAG
void BasicExpression::AddExpressionTreeToDAG(map<int,vector<int> >& DAG,
					     map<int,int>& VOp,
					     map<int,set<int> >& VColor,
					     map<int,int>& ExprSymm,
					     int eidx,
					     map<double,int,DoubleLessThan>&
					     ConstantSymm,
					     int varsymmclasses,
					     map<int,int>& VarSymm,
					     int maxcolorexpr,
					     int startnode,
					     int& topnode,
					     int totalops,
					     int level) {
  int thisnode = startnode + topnode;
  int op = GetOpType();
  int allconsts = ConstantSymm.size() + 1;
  int allvars = varsymmclasses;
  int allord = 2;
  int color = 0;
  int vi;

  if (IsConstant()) {
    // deal with constants
    color = computecolor(ExprSymm[eidx], level, op,
			 ConstantSymm[GetValue()], 0, 0,
			 maxcolorexpr, totalops, allconsts, allvars, allord);
    VColor[color].insert(thisnode);
  } else {

    // special cases: non-unit coefficients and exponents
    if (fabs(GetCoeff() - 1) > DAGTOLERANCE &&
	fabs(GetExponent() - 1) <= DAGTOLERANCE) {
      // if node is "coeff*node", transform it in three nodes
      // product node
      color = computecolor(ExprSymm[eidx], level, PRODUCT, 0, 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(thisnode);
      VOp[thisnode] = PRODUCT;
      // constant node
      level++;
      topnode++;
      color = computecolor(ExprSymm[eidx], level, CONST,
			   ConstantSymm[GetCoeff()], 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(startnode + topnode);
      DAG[thisnode].push_back(startnode + topnode);
      VOp[startnode + topnode] = CONST;
      // op node itself
      if (!IsVariable()) {
	topnode++;
	color = computecolor(ExprSymm[eidx], level, op, 0, 0, 0,
			     maxcolorexpr, totalops, allconsts, allvars,
			     allord);
	VColor[color].insert(startnode + topnode);
	DAG[thisnode].push_back(startnode + topnode);
	VOp[startnode + topnode] = op;
	// update thisnode
	thisnode = startnode + topnode;
      } else {
	vi = GetVarIndex();
	color = computecolor(0, 0, 0, 0, VarSymm[vi], 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(vi);
	DAG[thisnode].push_back(vi);
	VOp[vi] = VAR;
      }
    } else if (fabs(GetCoeff() - 1) <= DAGTOLERANCE &&
	       fabs(GetExponent() - 1) > DAGTOLERANCE) {
      // if node is "node^exponent", transform it in three nodes
      // power node
      color = computecolor(ExprSymm[eidx], level, POWER, 0, 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(thisnode);
      VOp[thisnode] = POWER;
      // constant node
      level++;
      topnode++;
      color = computecolor(ExprSymm[eidx], level, CONST,
			   ConstantSymm[GetExponent()], 0, 1, //1=second arg
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(startnode + topnode);
      DAG[thisnode].push_back(startnode + topnode);
      VOp[startnode + topnode] = CONST;
      // op node itself
      if (!IsVariable()) {
	topnode++;
	color = computecolor(ExprSymm[eidx], level, op, 0, 0, 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(startnode + topnode);
	DAG[thisnode].push_back(startnode + topnode);
	VOp[startnode + topnode] = op;
	// update thisnode
	thisnode = startnode + topnode;
      } else {
	vi = GetVarIndex();
	color = computecolor(0, 0, 0, 0, VarSymm[vi], 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(vi);
	DAG[thisnode].push_back(vi);
	VOp[vi] = VAR;
      }
    } else if (fabs(GetCoeff() - 1) > DAGTOLERANCE &&
	       fabs(GetExponent() - 1) > DAGTOLERANCE) {
      // if node is "coeff*(node^exponent)", transform it in 4 nodes
      // product node
      color = computecolor(ExprSymm[eidx], level, PRODUCT, 0, 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(thisnode);
      VOp[thisnode] = PRODUCT;
      // coefficient node
      level++;
      topnode++;
      color = computecolor(ExprSymm[eidx], level, CONST,
			   ConstantSymm[GetCoeff()], 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(startnode + topnode);
      DAG[thisnode].push_back(startnode + topnode);
      VOp[startnode + topnode] = CONST;
      // power node
      topnode++;
      color = computecolor(ExprSymm[eidx], level, POWER, 0, 0, 0,
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(startnode + topnode);
      DAG[thisnode].push_back(startnode + topnode);
      VOp[startnode + topnode] = POWER;
      thisnode = startnode + topnode;
      // exponent node
      level++;
      topnode++;
      color = computecolor(ExprSymm[eidx], level, CONST,
			   ConstantSymm[GetExponent()], 0, 1, // 1=second arg
			   maxcolorexpr, totalops, allconsts, allvars, allord);
      VColor[color].insert(startnode + topnode);
      DAG[thisnode].push_back(startnode + topnode);
      VOp[startnode + topnode] = CONST;
      // op node itself
      if (!IsVariable()) {
	topnode++;
	color = computecolor(ExprSymm[eidx], level, op, 0, 0, 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(startnode + topnode);
	DAG[thisnode].push_back(startnode + topnode);
	VOp[startnode + topnode] = op;
	// update thisnode
	thisnode = startnode + topnode;
      } else {
	vi = GetVarIndex();
	color = computecolor(0, 0, 0, 0, VarSymm[vi], 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(vi);
	DAG[thisnode].push_back(vi);
	VOp[vi] = VAR;
      }
    } else {
      if (IsVariable()) {
	vi = GetVarIndex();
	color = computecolor(0, 0, 0, 0, VarSymm[vi], 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VColor[color].insert(vi);
	VOp[vi] = VAR;
      } else {
	color = computecolor(ExprSymm[eidx], level, op, 0, 0, 0, maxcolorexpr,
			     totalops, allconsts, allvars, allord);
	VOp[thisnode] = op;
      }
    }

    // normal treatment of a node
    if (!IsVariable()) {
      VColor[color].insert(thisnode);
      int sz = GetSize();
      for(int i = 0; i < sz; i++) {
	if (GetNode(i)->IsVariable() &&
	    fabs(GetNode(i)->GetCoeff() - 1) <= DAGTOLERANCE &&
	    fabs(GetNode(i)->GetExponent() - 1) <= DAGTOLERANCE) {
	  // if child node is a "simple node" variable then deal with it here
	  // otherwise it's dealt with in the cases above
	  DAG[thisnode].push_back(GetNode(i)->GetVarIndex());
	} else {
	  topnode++;
	  DAG[thisnode].push_back(startnode + topnode);
	}
	GetNode(i)->AddExpressionTreeToDAG(DAG, VOp, VColor,
					   ExprSymm, eidx,
					   ConstantSymm,
					   allvars, VarSymm,
					   maxcolorexpr,
					   startnode, topnode,
					   totalops, level+1);
	// WORKING HERE: missing: order for noncommutative operators
      }
    }
  }

}

#ifdef JUSTCOPYITTOWHEREYOUNEEDIT
	cout << "///////// color assignment" << endl;
	cout << "  eidx = " << eidx << endl;
	cout << "  ExprSymm[" << eidx << "] = " << ExprSymm[eidx] << endl;
	cout << "  level = " << level << endl;
	cout << "  op = " << op << endl;
	cout << "  maxcolorexpr = " << maxcolorexpr << endl;
	cout << "  totalops = " << totalops << endl;
	cout << "  allconsts = " << allconsts << endl;
	cout << "  allvars = " << allvars << endl;
	cout << "  allord = " << allord << endl;
	cout << "  ---------" << endl;
	cout << "  thisnode = " << thisnode << endl;
	cout << "  op = " << op << endl;
	cout << "  color = " << color << endl;
	cout << "///////// end color assignment" << endl;
#endif

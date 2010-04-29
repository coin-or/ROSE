/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        operand.cxx                                            *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (operand class functions)   *
* History:     080928 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#include <string>
#if __GNUC__ >= 3
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

// constructors
Operand::Operand() : 
  oplabel(CONST), constant(0), varindex(NOVARIABLE), 
  coefficient(1), exponent(1), depconstant(NULL), depcoefficient(NULL), 
  depexponent(NULL), dependency(0), lb(-LARGE), ub(LARGE) { }

Operand::Operand(double t) : 
  oplabel(CONST), constant(t), varindex(NOVARIABLE),
  coefficient(1), exponent(1), depconstant(NULL), depcoefficient(NULL), 
  depexponent(NULL), dependency(0), lb(-LARGE), ub(LARGE) { }

Operand::Operand(const Int t) : 
  oplabel(CONST), constant(t), varindex(NOVARIABLE),
  coefficient(1), exponent(1), depconstant(NULL), depcoefficient(NULL), 
  depexponent(NULL), dependency(0), lb(-LARGE), ub(LARGE) { }

Operand::Operand(const Int t, const bool isvar) :   
  coefficient(1), exponent(1), depconstant(NULL), depcoefficient(NULL), 
  depexponent(NULL), dependency(0), lb(-LARGE), ub(LARGE) { 
  if (isvar) {
    // make it a variable
    oplabel = VAR;
    constant = 0;
    varindex = t;
  } else {
    // make it an operator label
    oplabel = (int) t;
    constant = 0;
    varindex = NOVARIABLE;
  }
}

// create an (empty) operator or a variable leaf and set coefficient
Operand::Operand(const double c, const Int t, string vn) : 
  coefficient(c), exponent(1), lb(-LARGE), ub(LARGE) { 
  // make it a variable
  oplabel = VAR;
  constant = 0;
  varindex = t;
  varname = vn;
  dependency = 0;
  depconstant = NULL;
  depcoefficient = NULL;
  depexponent = NULL;
}

// copy constructor
Operand::Operand(const Operand& t) : 
  oplabel(t.oplabel), constant(t.constant), depconstant(t.depconstant), 
  varindex(t.varindex), varname(t.varname), coefficient(t.coefficient), 
  depcoefficient(t.depcoefficient), exponent(t.exponent), 
  depexponent(t.depexponent), dependency(t.dependency),
  lb(t.lb), ub(t.ub) { }

// destructor
Operand::~Operand() { }

// Operand class methods:

string Operand::ToString(void) const {
#if __GNUC__ >= 3
  stringstream outbuf;  
#else
  strstream outbuf;
#endif
  string vn;
  if (GetCoeff() == 0) {
    // coefficient is 0
    outbuf << 0;
  } else if (GetOpType() == CONST) {
    // constant
    outbuf << GetValue();
  } else if (GetOpType() == VAR) {
    // variable
    if (GetCoeff() == 1) {
      int vi = GetVarIndex();
      if (vi == NOVARIABLE) {
	if (GetExponent() == 1) {
	  outbuf << NOTVARNAME;
	} else {
	  outbuf << NOTVARNAME << "^" << GetExponent();
	}
      } else {
	vn = GetVarName();
	if (GetExponent() == 1) {
	  //outbuf << vn << VNAMEIDXCHAR << vi;
	  outbuf << vn;
	} else {
	  //outbuf << vn << VNAMEIDXCHAR << vi << "^" << GetExponent();
	  outbuf << vn << "^" << GetExponent();
	}
      }
    } else {
      int vi = GetVarIndex();
      if (vi == NOVARIABLE) {
	if (GetExponent() == 1) {
	  outbuf << GetCoeff() << "*_" << NOTVARNAME;
	} else {
	  outbuf << GetCoeff() << "*" << NOTVARNAME << "^" 
                 << GetExponent();
	}
      } else {
	vn = GetVarName();
	if (GetExponent() == 1) {
	  //outbuf << GetCoeff() << "*" << vn << VNAMEIDXCHAR << vi;
	  outbuf << GetCoeff() << "*" << vn;
	} else {
	  //outbuf << GetCoeff() << "*" << vn << VNAMEIDXCHAR << vi << "^" 
	  //<< GetExponent();
	  outbuf << GetCoeff() << "*" << vn << "^" << GetExponent();
	}
      }
    }
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
double Operand::GetValue(void) const { 
  double ret = 0;
  if (oplabel == CONST && dependency == 0) {
    if (exponent == 1) 
      ret = coefficient * constant;
    else if (exponent == 2)
      ret = coefficient * constant * constant;
    else
      ret = coefficient * pow(constant, exponent);
  } else if (oplabel == CONST && dependency == 1 && depconstant) {
    ret = *depconstant;
  } else {
    ret = constant; 
  }
  return ret;
}

// just get the value in any case
double Operand::GetSimpleValue(void) const { 
  if (dependency == 1 && depconstant) {
    return *depconstant;
  } else {
    return constant; 
  }
}

// node bounds
double Operand::GetLB(void) const {
  return lb;
}
double Operand::GetUB(void) const {
  return ub;
}
void Operand::SetLB(double theLB) {
  lb = theLB;
}
void Operand::SetUB(double theUB) {
  ub = theUB;
}

// get variable index
Int Operand::GetVarIndex(void) const { return varindex; }

// get variable name
string Operand::GetVarName(void) const { return varname; }

// set operator type
void Operand::SetOpType(const int t) { oplabel = t; }
  
// set constant value
void Operand::SetValue(const double t) { oplabel = CONST; constant = t; }

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

void Operand::SetCoeff(const double coeff) {
  coefficient = coeff;
}

double Operand::GetCoeff(void) const {
  if (dependency == 2 && depcoefficient) {
    return *depcoefficient;
  } else {
    return coefficient;
  }
}

void Operand::SetExponent(const double expon) {
  exponent = expon;
}

double Operand::GetExponent(void) const {
  if (dependency == 3 && depexponent) {
    return *depexponent;
  } else {
    return exponent;
  }
}

void Operand::SetDependencyOnOperand(const int whichconstant, 
				     double** depvalue) {
  dependency = whichconstant + 1;
  switch(dependency) {
  case 1:
    depconstant = *depvalue;
    break;
  case 2:
    depcoefficient = *depvalue;
    break;
  case 3:
    depexponent = *depvalue;
    break;
  }
}

void Operand::EnforceDependencyOnOperand(void) {
  switch(dependency) {
  case 1:
    constant = *depconstant;
    break;
  case 2:
    coefficient = *depcoefficient;
    break;
  case 3:
    exponent = *depexponent;
    break;
  }
}  

void Operand::ConsolidateValue(void) {
  SetValue(GetValue());
  SetCoeff(1.0);
  SetExponent(1.0);
}

// is operand a zero constant?
bool Operand::IsZero(void) const {
  if (GetOpType() == CONST) {
    double c = GetValue();
    if (fabs(c) < Ev3NearZero()) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
  
// is operand a constant having value v?
bool Operand::HasValue(double v) const {
  if (GetOpType() == CONST) {
    double c = GetValue();
    double t1 = v + Ev3NearZero();
    double t2 = v - Ev3NearZero();
    if (c < t1 && c > t2) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}  

// is operand a negative constant?
bool Operand::IsLessThan(double v) const {
  if (GetOpType() == CONST) {
    if (GetValue() < v + Ev3NearZero()) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

// is operand a negative constant?
bool Operand::IsGreaterThan(double v) const {
  if (GetOpType() == CONST) {
    if (GetValue() > v - Ev3NearZero()) {
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
void Operand::SubstituteVariableWithConstant(int varindex, double c) {
  if (GetOpType() == VAR && GetVarIndex() == varindex) {  
    SetOpType(CONST);
    SetVarIndex(NOVARIABLE);
    double t;
    t = GetCoeff() * pow(c, GetExponent());
    SetCoeff(1);
    SetExponent(1);
    SetValue(t);    
  }
}


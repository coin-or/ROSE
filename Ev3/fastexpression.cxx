/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        fastexpression.cxx                                     *
* Source:      GNU C++                                                *
* Purpose:     simplified non-C++ expression trees for fast eval/diff *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
* History:     050909 0.0 work started                                *
***********************************************************************/

#include <iostream>
#include <map>
#include <cmath>
#include "common.h"
#include "fastexpression.h"

void DeleteFastEvalTree(FastEvalTree* fet) {
  if (fet) {
    DeleteFastEvalTreeRecursive(fet);
    delete fet;
  }
}

void DeleteFastEvalTreeRecursive(FastEvalTree* fet) {
  if (fet) {
    for(Int i = 0; i < fet->nodesize; i++) {
      DeleteFastEvalTreeRecursive(&fet->nodes[i]);
    }
    delete [] fet->nodes;
  }
}

double FastEval(FastEvalTree* fet, double* varvalues, Int vsize) {
  return FastEvalRecursive(fet, varvalues, NULL, vsize);
}

double FastEval(FastEvalTree* fet, double* varvalues, 
		std::map<int,int>& varmap, Int vsize) {
  return FastEvalRecursive(fet, varvalues, &varmap, vsize);
}

double FastEvalRecursive(FastEvalTree* fet, double* varvalues, 
			 std::map<int,int>* varmap, Int vsize) {
  using namespace std;
  double thecoeff = fet->depcoeff ? *(fet->depcoeff) : fet->coeff;
  double theexpon = fet->depexponent ? *(fet->depexponent) : fet->exponent;
  double thevalue = fet->depvalue ? *(fet->depvalue) : fet->value;
  int localindex;
  // evaluate expression
  if (fet->optype == CONST) {
    return thevalue;
  } else if (fet->optype == VAR) {
    if (varmap) {
      localindex = (*varmap)[fet->varindex];
    } else {
      localindex = fet->varindex;
    }
    assert(localindex <= vsize);
    assert(localindex > 0);
    double ret = 0;
    ret = varvalues[localindex - 1];
    ret = thecoeff * pow(ret, theexpon);
    return ret;
  } else { 
    Int i;
    double ret = 0;
    double tmp = 0;
    int op = fet->optype;
    assert(fet->nodesize > 0);
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
	cerr << "FastEvalRecursive: in [GetSize()!=2]: "
	     << "fractions should have just 2 operands, but going ahead "
	     << "anyway...\n";
      }
      for(i = 0; i < fet->nodesize; i++) {
	tmp = FastEvalRecursive(&fet->nodes[i], varvalues, varmap, vsize);
	if (i > 0 && tmp == 0) {
	  cerr << "FastEvalRecursive: division by zero not allowed (node "
	       << i << " evaluates to 0), setting to a large value" << endl;
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
	cerr << "FastEvalRecursive: in [GetSize()!=2]: "
	     << "powers should have just 2 operands, but going ahead "
	     << "anyway...\n";
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
	cerr << "FastEvalRecursive: log of negative not allowed ("
	     << "argument evaluates to < 0), taking absolute value" << endl;
	tmp = -tmp;
      } else if (tmp == 0) {
	cerr << "FastEvalRecursive: log of zero not allowed ("
	     << "argument evaluates to 0), setting to large negative value" 
	     << endl;
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
	cerr << "FastEvalRecursive: sqrt of negative not allowed, "
	     << "taking absolute value" << endl;
	tmp = -tmp;
      }
      ret = thecoeff * sqrt(tmp);
      break;
    }
    return ret;
  }
}

FastEvalTree* Diff(FastEvalTree* fet, Int varindex) {
  using namespace std;
  FastEvalTree* ret = NULL;
  cerr << "FastEvalTree|Diff not implemented\n";
  return ret;
}

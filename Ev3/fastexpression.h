/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        fastexpression.h                                       *
* Source:      GNU C++                                                *
* Purpose:     non-C++ expression trees for fast eval/diff            *
* History:     050909 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3FASTEXPRESSIONH__
#define __EV3FASTEXPRESSIONH__

#define RCS13 "$Id: fastexpression.h,v 1.2 2006/07/30 05:36:45 liberti Exp liberti $"

#include "common.h"

struct fevaltree {
  int optype;
  Int varindex;
  double coeff;
  double* depcoeff;
  double exponent;
  double* depexponent;
  double value;
  double* depvalue;
  struct fevaltree* nodes;
  Int nodesize;
};
typedef struct fevaltree FastEvalTree;

// derivative of fet w.r.t. variable index varindex in (separate) tree
FastEvalTree* Diff(FastEvalTree* fet, int varindex);
// deletion (use nonrecursive version)
void DeleteFastEvalTree(FastEvalTree* fet);
void DeleteFastEvalTreeRecursive(FastEvalTree* fet);
// Fast evaluators (for repeated evaluations -- use nonrecursive version)
double FastEval(FastEvalTree* fet, double* varvalues, Int vsize);
double FastEval(FastEvalTree* fet, double* varvalues, 
		std::map<int,int>& varmap, Int vsize);  
double FastEvalRecursive(FastEvalTree* fet, double* varvalues, 
			 std::map<int,int>* varmap, Int vsize);


#endif

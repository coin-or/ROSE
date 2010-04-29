/*
** Name:       solver.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    "Common" problem solver implementation file
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050922 work started
*/

#include <iostream>
#include <vector>
#include "problem.h"
#include "solver.h"

UserCut::UserCut(Expression e, double theLB, double theUB) : 
  Function(e), LB(theLB), UB(theUB) { 
  NonlinearPart.SetTo(Function->GetNonlinearPart());
  FunctionFET = Function->GetFastEvalTree();
  NonlinearPartFET = Function->GetFastEvalTree();
}

UserLinearCut::UserLinearCut() : 
  LB(-ROSEINFINITY), UB(ROSEINFINITY), Nonzeroes(0), memcount(NULL),
  VarIndices(NULL), Coeffs(NULL) { }

UserLinearCut::UserLinearCut(std::vector<std::pair<int, double> >& lininfo, 
			     double theLB, double theUB) :
  LB(theLB), UB(theUB) {
  using namespace std;
  int nzulcut = 0;
  for(vector<pair<int,double> >::iterator vpi = lininfo.begin();
      vpi != lininfo.end(); vpi++) {
    if (std::fabs(vpi->second) > EPSILONTOLERANCE) {
      nzulcut++;
    }
  }
  Nonzeroes = nzulcut;
  memcount = new int (1);
  VarIndices = new int [nzulcut + 1];
  Coeffs = new double [nzulcut + 1];
  VarIndices[0] = 0;
  Coeffs[0] = 0;
  int i = 1;
  for(vector<pair<int,double> >::iterator vpi = lininfo.begin();
      vpi != lininfo.end(); vpi++) {
    if (std::fabs(vpi->second) > EPSILONTOLERANCE) {
      VarIndices[i] = vpi->first;
      Coeffs[i] = vpi->second;
      i++;
    }
  }
}

UserLinearCut::UserLinearCut(int* varindices, double* coeffs, int asize, 
			     double theLB, double theUB) :
  LB(theLB), UB(theUB) {

  using namespace std;
  int nzulcut = 0;
  for(int i = 0; i < asize; i++) {
    if (std::fabs(coeffs[i]) > EPSILONTOLERANCE) {
      nzulcut++; 
    }
  }
  Nonzeroes = nzulcut;
  memcount = new int (1);
  VarIndices = new int [nzulcut + 1];
  Coeffs = new double [nzulcut + 1];
  VarIndices[0] = 0;
  Coeffs[0] = 0;
  int k = 1;
  for(int i = 0; i < asize; i++) {
    if (std::fabs(coeffs[i]) > EPSILONTOLERANCE) {
      VarIndices[k] = varindices[i];
      Coeffs[k] = coeffs[i];
      k++;
    }
  }
}

UserLinearCut::UserLinearCut(const UserLinearCut& uc) :
  LB(uc.LB), UB(uc.UB), Nonzeroes(uc.Nonzeroes), VarIndices(uc.VarIndices),
  Coeffs(uc.Coeffs), memcount(uc.memcount) {
  if (memcount) {
    (*memcount)++;
  }
}

UserLinearCut& UserLinearCut::operator=(const UserLinearCut& uc) {
  LB = uc.LB;
  UB = uc.UB;
  Nonzeroes = uc.Nonzeroes;
  VarIndices = uc.VarIndices;
  Coeffs = uc.Coeffs;
  memcount = uc.memcount;
  (*memcount)++;
  return *this;
}

UserLinearCut::~UserLinearCut(void) {
  if (memcount) {
    if (*memcount == 0) {
      if (VarIndices) {
	delete [] VarIndices;
      }
      if (Coeffs) {
	delete [] Coeffs;
      }
      delete memcount;
    } else {
      *memcount--;
    }
  }
}


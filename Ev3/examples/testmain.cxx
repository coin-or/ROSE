/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        testmain.cxx                                           *
* Source:      GNU C++                                                *
* Purpose:     main driver program for Expression v3                  *
* History:     010517 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#define RCS6 "$Id: testmain.cxx,v 1.2 2001/07/23 01:50:06 liberti Exp liberti $"

#include <iostream>

#include "expression.h"
#include "parser.h"

using namespace std;

int main(int argc, char** argv) {

  Expression e1(1.0);
  Expression e2(2.0);
  Expression ex; // true --> 1 is a variable index
  ex->SetOpType(VAR);
  ex->SetVarIndex(1);
  Expression ey;
  ex->SetOpType(VAR);
  ex->SetVarIndex(2);
  Expression ex2(0.5);
  Expression ef;
  ef->SetOpType(SUM);

  ef->AddNode(ex);
  ef->AddNode(ey);
  ef->AddNode(e2);

  ef->Debug();

  cout << "comparison of ef, ey: ";
  cout << (ef == ey) << endl;

  cout << "sum of e2, e1: ";
  e2 + e1;
  cout << e2->GetValue() << endl;

  cout << "sum of ex2, ex: ";
  ex2 + ex;
  cout << ex2->GetCoeff() << endl;

  cout << "sum of ef, e1: ";
  cout << ef->GetSize() << ", ";
  ef + e1;
  cout << ef->GetSize() << endl;

  cout << "sum of e1, ef: ";
  e1 + ef;
  cout << e1->GetSize() << endl;

  cout << "comparison of e1, e1: ";
  cout << (e1 == e1) << endl;

  return 0;
}

/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        test.cxx                                               *
* Source:      GNU C++                                                *
* Purpose:     main driver program for Expression v3                  *
* History:     021218 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#include "expression.h"
#include "parser.h"

#include <iostream>
#include <cstring>

using namespace std;

#define BUFSIZE 8192

int main(int argc, char** argv) {

  ExpressionParser p;
  Expression e;
  int nerr;

  string es = "-x1^2 -x2^2 - x3^2 - x1*x2";
  e = p.Parse(es.c_str(), nerr);
  Simplify(&e);
  cout << es << " simplifies to: " << endl;
  cout << "  " << e->ToString() << endl;
  cout << "derivatives: " << endl;
  cout << "  w.r.t. x1 = " << Diff(e, 1)->ToString() << endl;
  cout << "  w.r.t. x2 = " << Diff(e, 2)->ToString() << endl;
  cout << "  w.r.t. x3 = " << Diff(e, 3)->ToString() << endl;

  return 0;
}

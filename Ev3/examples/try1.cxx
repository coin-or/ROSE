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

  char buffer[BUFSIZE];

  ExpressionParser p;
  Expression e, f;
  int nerr;
  bool ischanged;

  e = p.Parse("x+y*z-2*x+x", nerr);
  f = p.Parse("z*y", nerr);

  Simplify(&e);
  Simplify(&f);

  cout << e->ToString() << endl;
  cout << f->ToString() << endl;
  cout << (e == f) << endl;

  return 0;
}

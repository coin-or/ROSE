/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        examplemain.cxx                                        *
* Source:      GNU C++                                                *
* Purpose:     main driver program for Expression v3                  *
* History:     010721 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#define RCS12 "$Id: examplemain.cxx,v 1.6 2002/12/09 18:49:31 liberti Exp liberti $"

#include <iostream>
#include "expression.h"
#include "parser.h"

using namespace std;

int main(int argc, char** argv) {

  // test parser
  ExpressionParser p;       // create the parser object
  p.SetVariableID("x", 1);  // configure it with variables x...
  p.SetVariableID("y", 2);  // ... and y
  p.SetVariableID("z", 3);
  int nerr = 0;             // nerr contains the number of parser errors
  Expression pe(p.Parse("log(2*x*y)+sin(z)", nerr)); // parse
  cout << "parsing errors: " << nerr << endl;
  cout << "f = " << pe->ToString() << endl;
  //cout << "tree is\n" << pe->PrintTree(0,2) << endl;
  double val[3] = {2, 3, 1};
  cout << "eval(2,3,1):   " << pe->Eval(val, 3) << endl;
  cout << "numeric check: " << ::log(2*2*3)+::sin(1) << endl;
  // test diff
  Expression de1 = Diff(pe, 1);
  cout << "df/dx = " << de1->ToString() << endl;
  Expression de2 = Diff(pe, 2);
  cout << "df/dy = " << de2->ToString() << endl;
  Expression de3 = Diff(pe, 3);
  cout << "df/dz = " << de3->ToString() << endl;
  // test TrigSimp
  Expression t(p.Parse("(sin(x+y+z))^2+(cos(x+y+z))^2", nerr));
  cout << t->ToString() << endl;
  TrigSimp(t);
  cout << "after TrigSimp: " << t->ToString() << endl;
  return 0;
}

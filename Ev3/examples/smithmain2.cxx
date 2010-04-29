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
#include <sys/time.h>
#include <cstdlib>

#define BUFSIZE 8192
#define DEFAULTVALUE 100

using namespace std;

int main(int argc, char** argv) {

  ExpressionParser p;
  Expression e, f;
  int nerr;

  // e = p.Parse("(x-2)^2 + (y-3)^2 + (z+1)^2", nerr);
  // e=p.Parse("sin (13*(x-y)) * (x*x + y*exp (x)) + 3 * log (1 + x^4)", nerr);
  // e = p.Parse("(x1+x3)*x2 - 2*x1*(x3-1) +5*x1 + 3*x2 + x3^2", nerr);
  // e = p.Parse("(x+y)*(x-z)*x*z^2", nerr);
  e = p.Parse("x*y", nerr);
  Simplify(&e);
  cout << e->ToString() << endl;

  // test smith's std form
  /*
  int vi = 4;
  string vn = "w";
  vector<int> oplabels;
  oplabels.push_back(SUM);
  oplabels.push_back(PRODUCT);
  vector<Expression> schemata;
  Expression x(1,1,"x");
  x->SetExponent(2);
  schemata.push_back(x);
  vector<Expression> defcons;
  e->SmithStandardForm(vi, vn, oplabels, schemata, defcons);
  for(int i = 0; i < defcons.size(); i++) {
    cout << vn << vi + i << " = " << defcons[i]->ToString() << endl;
  }
  cout << e->ToString() << endl;
  */

  // test replacement
  Expression needle1(1.0,1,"x");
  Expression replace1(0.0);
  for(int i = 0; i < 4; i++) {
    Expression xi(1.0,i+3, "x");
    replace1 = replace1 + xi;
  }
  Expression needle2(1.0,2,"y");
  Expression replace2(0.0);
  for(int i = 0; i < 4; i++) {
    Expression yi(1.0,i+7, "y");
    replace2 = replace2 + yi;
  }
  e->ReplaceSubexpression(needle1, replace1);
  e->ReplaceSubexpression(needle2, replace2);
  cout << e->ToString() << endl;
  while(e->DistributeProductsOverSums());
  cout << e->ToString() << endl;
  int vi = 11;
  string vn = "w";
  vector<int> oplabels;
  oplabels.push_back(SUM);
  oplabels.push_back(PRODUCT);
  vector<Expression> schemata;
  Expression x(1,1,"x");
  x->SetExponent(2);
  schemata.push_back(x);
  vector<Expression> defcons;
  e->SmithStandardForm(vi, vn, oplabels, schemata, defcons);
  for(int i = 0; i < defcons.size(); i++) {
    cout << vn << "_" << vi + i << " = " << defcons[i]->ToString() << endl;
  }
  cout << e->ToString() << endl;
  return 0;
}

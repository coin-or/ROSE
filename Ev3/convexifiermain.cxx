/**********************************************************************
* Author:      Leo Liberti & Sonia Cafieri                                           *
* Name:        test.cxx                                               *
* Source:      GNU C++                                                *
* Purpose:     main driver program for Expression v3                  *
* History:     080603 work started                                *
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

int main(int argc, char** argv) {
  using namespace std;
  ExpressionParser p;
  Expression e, f;
  int nerr;


  e = p.Parse("x1*x2", nerr);
  Simplify(&e);
  cout << e->ToString() << endl;

  int vi = 1;
  string vn = "w";
  
  Expression x(1, 1, "x"), y(1, 2, "y");
  Expression schema = x * y;

/*

  // replace by product
  f = e->ReplaceByOperator(vi, vn, PRODUCT);
  while(!f->IsZero()) {
cout <<"replace by product " << endl;
    cout << vn << "[" << vi << "] = " << f->ToString() << endl;
    vi++;
    f = e->ReplaceByOperator(vi, vn, PRODUCT);
  }

  // replace by power operator
  f = e->ReplaceByOperator(vi, vn, POWER);
  while(!f->IsZero()) {
cout <<"replace by power " << endl;
    cout << vn << "[" << vi << "] = " << f->ToString() << endl;
    vi++;
    f = e->ReplaceByOperator(vi, vn, POWER);
  }

*/

cout << "replace by schema" << endl;
    f = e->ReplaceBySchema(vi, vn, schema);

  cout << e->ToString() << endl;
  cout << f->ToString() << endl;

    double vlb[100], vub[100]; 
    vlb[0] = 1.; vub[0] = 2.;
    vlb[1] = 1.; vub[1] = 2.;
 
    // xj -> x, xL1, xU1
    // xk -> y, xL2, xU2
    //  Expression x(1.0, x1, "x");
    //  Expression y(1.0, x2, "y");

      Expression xL1(vlb[0]);
      Expression xU1(vub[0]);
      Expression xL2(vlb[1]);
      Expression xU2(vub[1]);
      Expression c1 = xL1*y + xL2*x - xL1*xL2;
      Expression c2 = xU1*y + xU2*x - xU1*xU2;
      Expression c3 = xL1*y + xU2*x - xL1*xU2;
      Expression c4 = xU1*y + xL2*x - xU1*xL2;
      Simplify(&c1);
      Simplify(&c2);
      Simplify(&c3);
      Simplify(&c4);

      // add constraint
      vn = "c";
    //  AddConstraint(vn, c1, 0.0, 0.0);

  return 0;
}

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

int main(int argc, char** argv) {

  ExpressionParser p;
  Expression e, f;
  int nerr;


  // e = p.Parse("(x-2)^2 + (y-3)^2 + (z+1)^2", nerr);  
  // e=p.Parse("sin (13*(x-y)) * (x*x + y*exp (x)) + 3 * log (1 + x^4)", nerr);
  // e = p.Parse("(x1+x3)*x2 - 2*x1*(x3-1) +5*x1 + 3*x2 + x3^2", nerr);
  //e = p.Parse("(x+y)*(x-z)*x*z^2", nerr);
  e = p.Parse("x*y", nerr);
  Simplify(&e);
  cout << e->ToString() << endl;

  //int vi = 3;
  int vi = 4;
  string vn = "w";
  
  //Expression x(1, 1, "x1"), y(1, 2, "x2");
  Expression schema = x1 * y2;
    f = e->ReplaceBySchema(vi, vn, schema);
  cout << e->ToString() << endl;
  cout << f->ToString() << endl;

/*
  // replace by sum
  f = e->ReplaceByOperator(vi, vn, SUM);
  while(!f->IsZero()) {
cout <<"replace by sum " << "vi=" << vi <<endl;
    cout << vn << "[" << vi << "] = " << f->ToString() << endl;
    vi++;
    f = e->ReplaceByOperator(vi, vn, SUM);
  }
  vi++;

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

  // replace by operand
cout << "replace by schema" << endl;
  Expression x(1,1,"x");
  x->SetExponent(2);
  // replace by power operator
  f = e->ReplaceBySchema(vi, vn, x);
  while(!f->IsZero()) {
cout <<"replace by operator " << endl;
    cout << vn << "[" << vi << "] = " << f->ToString() << endl;
    vi++;
    f = e->ReplaceBySchema(vi, vn, x);
  }
  
  cout << e->ToString() << endl;
  
*/

  return 0;
}

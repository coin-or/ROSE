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


  e = p.Parse("x1*x2*x3*x4+x5*x6*x7*x8", nerr);
  Simplify(&e);
  cout << e->ToString() << endl;

  int vi = 1;
  string vn = "w";
  
  Expression x(1, 1, "x"), y(1, 2, "y"), z(1, 3, "z"),t(1, 4, "t");
  Expression schema = x * y *z *t;


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



    f = e->ReplaceBySchema(vi, vn, schema);
 // cout << e->ToString() << endl;
 // cout << f->ToString() << endl;
cout << "replace by schema" << endl;

  //Expression x(1,1,"x");
  //x->SetExponent(2);
  //f = e->ReplaceBySchema(vi, vn, x);
/*
  while(!f->IsZero()) {
    cout << vn << "[" << vi << "] = " << f->ToString() << endl;
    vi++;
    f = e->ReplaceBySchema(vi, vn, x);
  }
*/  
  cout << e->ToString() << endl;
  cout << f->ToString() << endl;
  


  return 0;
}

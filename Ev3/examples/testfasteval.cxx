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

using namespace std;

#define BUFSIZE 8192
#define DEFAULTVALUE 100

int main(int argc, char** argv) {

  ExpressionParser p;
  Expression e, f;
  int nerr;

  //  e = p.Parse("(x-2)^2 + (y-3)^2 + (z+1)^2", nerr);
  e = p.Parse("sin (13*(x-y)) * (x*x + y*exp (x)) + 3 * log (1 + x^4)", nerr);
  Simplify(&e);
  cout << e->ToString() << endl;

  int n = DEFAULTVALUE;
  if (argc > 1) {
    n = atoi(argv[1]);
    if (n == 0) {
      n = DEFAULTVALUE;
    }
  }

  // initialize randomizer
  struct timeval theTV;
  struct timezone theTZ;
  gettimeofday(&theTV, &theTZ);
  srand48(theTV.tv_usec);

  double x[3];
  double t;
  cout << "x\t\ty\t\tz\t\tf\tfeval\n";

  for(int i = n; i > 0; i--) {
    x[0] = (double) drand48();
    x[1] = (double) drand48();
    x[2] = (double) drand48();
    t = e->Eval(x, 2);
    cout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << t;
    t = e->FastEval(x, 2);
    cout << "\t" << t << endl;
  }

  return 0;
}

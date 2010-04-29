/**********************************************************************
* Author:      Franco Raimondi                                        *
* Name:        mbf.cxx                                                *
* Source:      GNU C++                                                *
* Purpose:     Test program for bf.cxx                                *
* History:     040609 work started                                    *
***********************************************************************/

#include "bf.h"

#include <iostream>
#include <cstring>

#define BUFSIZE 8192

int main(int argc, char** argv) {

  BasicBoolForm *a,*b,*c,*d;
  string s1 = "v1";
  string s2 = "v2";
  string s3 = "v3";
  a = new BasicBoolForm(1,s1);
  b = new BasicBoolForm(2,s2);
  c = new BasicBoolForm(3,s3);
  d = new BasicBoolForm(false);

  BoolForm e,f;
  e = !( (*a) + !(*b) ) * !(*c);
  f = !e;

  cout << "  " << e->ToString() << endl;
  cout << "  " << f->ToString() << endl;
  f = PushNeg(f);

  cout << "  " << f->ToString() << endl;
  return 0;
}
